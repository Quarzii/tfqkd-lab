"""Focused behavior checks for source fidelity and finite-band boundaries."""

import copy
import unittest
import numpy as np

from tfqkd.config import load
from tfqkd.integration import PhaseIntegral, frequency_grid
from tfqkd.reference import author_functions
from tfqkd.spectra import components, detection


class StageOneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c = load()
        cls.f = frequency_grid(cls.c, cls.c["grid"]["author_curve_points"])
        cls.ref = author_functions()

    def test_all_scenario_spectra_against_authors(self):
        for scenario in self.c["scenarios"]:
            with self.subTest(scenario=scenario["name"]):
                ref, _, _ = self.ref["calc_spectra"](scenario["common"], scenario["cavity"], scenario["stabilized"])
                expected = ref(self.f, self.c["operation"]["LB_km"], scenario["delta_L_km"])
                np.testing.assert_allclose(components(self.f, scenario, self.c)["total"], expected,
                                           rtol=self.c["validation"]["reference_rtol"], atol=0)

    def test_detection_once_and_K_only_common_fiber(self):
        changed = copy.deepcopy(self.c)
        changed["physics"]["K"] = 2
        for scenario in self.c["scenarios"]:
            original = components(self.f, scenario, self.c)
            revised = components(self.f, scenario, changed)
            np.testing.assert_array_equal(original["detection"], revised["detection"])
            np.testing.assert_array_equal(original["laser"], revised["laser"])
            # Eq. 5 versus Eq. 7, bertaina2024: changing K affects only common fiber.
            expected_fiber = original["fiber"] / 2 if scenario["common"] else original["fiber"]
            np.testing.assert_array_equal(revised["fiber"], expected_fiber)
            if scenario["stabilized"]:
                np.testing.assert_array_equal(original["detection"], detection(self.f, self.c["physics"]))

    def test_integral_endpoint_and_partial_bin(self):
        # Eq. 4, bertaina2024: constant synthetic PSD gives area of the interval.
        integral = PhaseIntegral([1., 2., 4.], [2., 2., 2.])
        np.testing.assert_allclose(integral.above([1., 1.5, 2., 3., 4., 5.]), [6., 5., 4., 2., 0., 0.])
        with self.assertRaises(ValueError):
            integral.above(0.5)

    def test_threshold_and_cap(self):
        # Eq. 4 and Sec. V, bertaina2024: analytic constant-PSD threshold example.
        integral = PhaseIntegral([1., 2., 4.], [2., 2., 2.])
        tau, reason = integral.operating_time(2., 1.)
        self.assertAlmostEqual(tau, 0.5)
        self.assertEqual(reason, "threshold")
        tau, reason = integral.operating_time(2., 0.4)
        self.assertEqual(tau, 0.4)
        self.assertEqual(reason, "capped")

    def test_bad_input_rejected(self):
        for f, psd in [([1., 1.], [1., 1.]), ([1., 2.], [1., -1.]), ([1., 2.], [1., np.nan])]:
            with self.assertRaises(ValueError):
                PhaseIntegral(f, psd)


if __name__ == "__main__":
    unittest.main()
