"""Phase-only physical checks; numerical tolerances are not apparatus coefficients."""
import unittest
import numpy as np
from scipy.special import sici

from tfqkd.mp_phase import (
    IntegrationGrid, components_from_config, component_variance,
    phase_difference_variance, tracking_band_sensitivity,
    fiber_component, laser_components, fiber_short_interval,
    equivalent_fiber_drift, zhang_x_error, zhou_x_error, discrete_average,
)
from tfqkd.spectra import free_laser, stabilized_laser


class ModePairingPhaseTests(unittest.TestCase):
    def setUp(self):
        # Numerical test grid; physical inputs below are explicit source examples.
        self.grid = IntegrationGrid(1e-8, 1e8, 16385, 16)
        self.fiber = {"length_km": 201.86, "l": 44., "fc1_hz": 100.}  # Table S1, Zhou; Table III, Bertaina.
        self.laser = {"model": "free", "r3": 3e6, "r2": 3e2, "fc_hz": 2e6}  # Table III, Bertaina.

    def test_white_frequency_normalization(self):
        c = laser_components("A", {"model": "lorentzian_white_frequency", "linewidth_hz": 100})[0]
        t = np.array([1e-5, 1e-4, 1e-3])  # Numerical intervals, not physical defaults.
        d, lower = component_variance(c, t, self.grid)
        # Eq. 4, didomenico2010: the closed finite-band integral, with h0=linewidth/pi.
        si, _ = sici(2 * np.pi * self.grid.upper_hz * t)
        expected = 4 * (100 / np.pi) * (np.pi * t * si
                                        - np.sin(np.pi * self.grid.upper_hz * t)**2 / self.grid.upper_hz)
        np.testing.assert_allclose(d, expected, rtol=1e-4)
        self.assertEqual(lower, 0)

    def test_small_interval_and_negative_cubic_correction(self):
        c = fiber_component("fiber", self.fiber)
        t = np.array([1e-7, 1e-6, 1e-5])
        d, _ = component_variance(c, t, self.grid)
        leading = fiber_short_interval(t, 44, 201.86, 100)
        corrected = fiber_short_interval(t, 44, 201.86, 100, next_order=True)
        self.assertTrue(np.all(d < leading))
        self.assertTrue(np.all(np.abs(d - corrected) < np.abs(d - leading)))
        np.testing.assert_allclose(d[:2], corrected[:2], rtol=2e-6)

    def test_no_round_trip_factor_and_arm_partition(self):
        a = fiber_component("A", dict(self.fiber, length_km=100.93))  # Table S1, zhou2023_async.
        b = fiber_component("B", dict(self.fiber, length_km=100.93))
        combined = fiber_component("single", self.fiber)
        t = np.array([1e-6, 1e-4])
        result = phase_difference_variance([a, b], t, self.grid)
        expected, _ = component_variance(combined, t, self.grid)
        np.testing.assert_allclose(result["variance_rad2"], expected, rtol=1e-14)
        np.testing.assert_allclose(result["components_rad2"]["A"], expected / 2, rtol=1e-14)

    def test_flicker_requires_explicit_window(self):
        comps = laser_components("A", self.laser)
        with self.assertRaisesRegex(ValueError, "diverges.*T_track"):
            phase_difference_variance(comps, 1e-5, self.grid)

    def test_only_divergent_terms_are_cut(self):
        comps = laser_components("A", self.laser) + [fiber_component("fiber", self.fiber)]
        result = tracking_band_sensitivity(comps, [1e-6, 1e-4], self.grid, .01)  # Explicit test band convention.
        self.assertEqual(result["nominal"]["cut_components"], ["A.r3"])
        self.assertEqual(result["nominal"]["lower_bounds_hz"], {"A.r3": 100., "A.r2": 0., "fiber": 0.})
        for name in ("A.r2", "fiber"):
            np.testing.assert_array_equal(result["shorter_x10"]["components_rad2"][name],
                                          result["longer_x10"]["components_rad2"][name])
        self.assertTrue(np.all(result["shorter_x10"]["variance_rad2"] < result["nominal"]["variance_rad2"]))
        self.assertTrue(np.all(result["nominal"]["variance_rad2"] < result["longer_x10"]["variance_rad2"]))

    def test_tracking_cutoff_below_numerical_floor_keeps_full_band(self):
        # User-approved cutoff convention; compare two quadratures over exactly the same band.
        component = laser_components("A", self.laser)[0]
        times = np.array([1e-6, 1e-4])
        fine_floor, _ = component_variance(component, times, self.grid, .01)
        coarse_floor_grid = IntegrationGrid(1000., 1e8, 16385, 16)
        covered, lower = component_variance(component, times, coarse_floor_grid, .01)
        self.assertEqual(lower, 100.)
        np.testing.assert_allclose(covered, fine_floor, rtol=2e-6)

    def test_laser_decomposition_matches_unchanged_tf_functions(self):
        f = self.grid.frequencies()
        np.testing.assert_allclose(sum(c.psd(f) for c in laser_components("A", self.laser)),
                                   free_laser(f, self.laser), rtol=1e-14)
        # Table III, bertaina2024; explicit reference parameters, no new presets.
        p = dict(self.laser, model="stabilized", C4=.5, C3=0., C2=2e-3,
                 B_hz=300e3, gamma=.1, delta=10.)
        components = laser_components("A", p)
        np.testing.assert_allclose(sum(c.psd(f) for c in components), stabilized_laser(f, p), rtol=1e-14)
        result = phase_difference_variance(components, 1e-5, self.grid, .01)
        self.assertEqual(result["cut_components"], ["A.C4"])

    def test_geometry_ignores_irrelevant_K(self):
        c = {"laser_A": self.laser, "laser_B": self.laser,
             "fiber_A": self.fiber, "fiber_B": self.fiber, "K": 4}
        original = phase_difference_variance(components_from_config(c), 1e-5, self.grid, .01)
        c["K"] = 2
        changed = phase_difference_variance(components_from_config(c), 1e-5, self.grid, .01)
        np.testing.assert_array_equal(original["variance_rad2"], changed["variance_rad2"])

    def test_known_raw_floor_and_frequency_offset(self):
        # Eq. B10, zhang2025, finite M=32; values from independent research check.
        self.assertAlmostEqual(float(zhang_x_error(0, 32)), .25160328721394853, places=14)
        # Eq. 2, zhou2023_async; V2=.46 and 1-kHz frequency from Fig. 3(a).
        e = zhou_x_error([0, 0], [0, .0005], .46, 1000)
        np.testing.assert_allclose(e, [.27, .73], rtol=1e-14)
        self.assertAlmostEqual(float(zhang_x_error(1e3, 32)), .5, places=14)

    def test_zero_interval_and_square_root_length(self):
        d, _ = component_variance(fiber_component("fiber", self.fiber), [0, 1e-6], self.grid)
        self.assertEqual(d[0], 0)
        # Eq. 6, bertaina2024, derived drift: doubling total length multiplies sigma by sqrt(2).
        ratio = equivalent_fiber_drift(44, 2 * 201.86, 100) / equivalent_fiber_drift(44, 201.86, 100)
        self.assertAlmostEqual(ratio, np.sqrt(2), places=14)

    def test_discrete_observed_count_average(self):
        # Table IV, zhu2023: actual four-bin QKD pair/error counts at 202 km.
        count = np.array([261923, 146520, 67458, 30799])
        errors = np.array([79234, 45185, 21410, 10150])
        # Eq. D5, zhang2025, discrete count-weighted mean; arithmetic using Table IV counts.
        expected = sum(errors) / sum(count)
        self.assertAlmostEqual(discrete_average(errors / count, count), expected, places=14)

    def test_input_rejections(self):
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            component_variance(fiber_component("fiber", self.fiber), -1, self.grid)
        with self.assertRaisesRegex(ValueError, "above"):
            component_variance(laser_components("A", self.laser)[0], 1e-5, self.grid, 1e-9)
        with self.assertRaisesRegex(ValueError, "weights"):
            discrete_average([.25, .3], [0, 0])
        with self.assertRaisesRegex(ValueError, "V2"):
            zhou_x_error(0, 1e-6, 1, 0)
        with self.assertRaises(KeyError):
            laser_components("A", {"model": "free", "r2": 300, "fc_hz": 2e6})


if __name__ == "__main__":
    unittest.main()
