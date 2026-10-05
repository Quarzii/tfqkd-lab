"""Independent normalization and inference checks of the new MP key module.

Synthetic photon yields and timing below are numerical fixtures, not apparatus
coefficients or defaults. Published inputs are separately checked by validator.
"""
import copy
import unittest

import numpy as np

from tfqkd.mp_key import (binary_entropy, pairing_rate, iid_pair_intervals,
                         phase_average, decoy_estimates, key_rate,
                         rate_per_second, calculate, markdown)


def photon_gains(intensities, error=False):
    """Source Poisson mixture with only 11,20,02 photons, known yields/errors."""
    x = np.asarray(intensities)[:, None]
    y = np.asarray(intensities)[None, :]
    # Eq.54/70, zeng2022 supplement: deliberately synthetic known photon yields.
    y11 = .02
    y20 = y02 = .02
    if error:
        # Eq.85, zeng2022 supplement: e0=1/2 for same-side photons; synthetic e11=.1.
        y11 *= .1
        y20 *= .5
        y02 *= .5
    return np.exp(-x-y) * (x*y*y11 + x*x/2*y20 + y*y/2*y02)


def configuration():
    # Table I/III and Eq.3, zhang2025: explicit mu/nu and linewidth/drift example.
    mu, nu = .234, .020
    return {
        'protocol': dict(mu=mu, nu=nu, click_probability=.001,
                         max_pairing_gap=128, z_pair_fraction=.1, e_z=.001, f_ec=1.06),
        'timing': dict(active_round_clock_hz=1.25e9, reference_fraction=1/8, recovery_fraction=0),
        'phase': dict(mode='published_drift', sigma_L_rad_s=4000., linewidth_hz=100.,
                      phase_slices=32, intrinsic_x_error=0.),
        'intervals': dict(mode='iid_zeng'),
        'decoy': dict(Q_Z=photon_gains([0, nu, mu]).tolist(),
                      Q_X=photon_gains([0, 2*nu, 2*mu]).tolist(),
                      QE_X=photon_gains([0, 2*nu, 2*mu], error=True).tolist(),
                      error_mode='observed'),
    }


class ModePairingKeyTests(unittest.TestCase):
    def test_entropy_endpoints_and_symmetry(self):
        np.testing.assert_allclose(binary_entropy([0, .5, 1]), [0, 1, 0], atol=0)
        self.assertAlmostEqual(float(binary_entropy(.13)), float(binary_entropy(.87)), places=14)

    def test_pairing_closed_limits(self):
        # Eq.4, zeng2022: exact single-gap and infinite-gap limits.
        self.assertEqual(pairing_rate(0, 1), 0)
        self.assertEqual(pairing_rate(1, 1), .5)
        p = .13
        self.assertAlmostEqual(pairing_rate(p, 1), p*p/(1+p), places=15)
        self.assertAlmostEqual(pairing_rate(p, 100000), p/2, places=15)

    def test_pairing_tiny_probability_is_resolved(self):
        # Eq.4, zeng2022, pL << 1: r_p ~ L*p^2, numerical stress only.
        p, gap = 1e-20, 100
        self.assertAlmostEqual(pairing_rate(p, gap)/(gap*p*p), 1., places=14)

    def test_geometric_gaps_and_endpoint(self):
        t, w = iid_pair_intervals(.1, 3, 10)
        # Eq.15/16, zeng2022: independently evaluated conditional probabilities.
        np.testing.assert_allclose(w, np.array([1, .9, .81])/2.71, rtol=1e-14)
        np.testing.assert_allclose(t, [.1, .2, .3], rtol=1e-14)
        _, w = iid_pair_intervals(1, 3, 10)
        np.testing.assert_array_equal(w, [1, 0, 0])

    def test_decoy_reconstructs_known_photon_error_without_floor_subtraction(self):
        c = configuration()
        p, d = c['protocol'], c['decoy']
        e = decoy_estimates(p['mu'], p['nu'], d['Q_Z'], d['Q_X'], d['QE_X'])
        self.assertAlmostEqual(e['Y_Z_11_lower_raw'], .02, places=13)
        self.assertAlmostEqual(e['Y_X_11_lower_raw'], .02, places=13)
        self.assertAlmostEqual(e['q_11_lower'], .5, places=13)
        self.assertAlmostEqual(e['e_11_x_upper'], .1, places=13)
        raw = d['QE_X'][1][1]/d['Q_X'][1][1]
        self.assertAlmostEqual(raw, .3, places=14)
        self.assertNotAlmostEqual(raw-.25, e['e_11_x_upper'], places=8)

    def test_decoy_common_normalization_invariance(self):
        c = configuration()
        p, d = c['protocol'], c['decoy']
        arguments = [np.asarray(d[k]) for k in ('Q_Z', 'Q_X', 'QE_X')]
        base = decoy_estimates(p['mu'], p['nu'], *arguments)
        scaled = decoy_estimates(p['mu'], p['nu'], *(a*17 for a in arguments))
        for key in ('q_11_lower', 'e_11_x_upper'):
            self.assertAlmostEqual(base[key], scaled[key], places=13)

    def test_unbounded_privacy_does_not_recover_key_above_half(self):
        c = configuration()
        p, d = c['protocol'], c['decoy']
        d['QE_X'] = (np.asarray(d['Q_X'])*.9).tolist()
        estimate = decoy_estimates(p['mu'], p['nu'], d['Q_Z'], d['Q_X'], d['QE_X'])
        self.assertAlmostEqual(estimate['e_11_x_upper'], .9, places=12)
        self.assertEqual(estimate['e_11_x_privacy_worst_case'], .5)
        self.assertTrue(estimate['notes'])
        self.assertEqual(calculate(c)['key_bps'], 0)

    def test_count_to_round_to_pair_to_second(self):
        # Eq.7, zeng2022: perfect fixture with exactly r_s*r_p bits per round.
        r = key_rate(.1, 3, .2, 1., 0., 0., 1.)
        self.assertAlmostEqual(r['bits_per_quantum_round'], r['r_p']*.2, places=15)
        self.assertAlmostEqual(r['bits_per_potential_pair'], 2*r['bits_per_quantum_round'], places=15)
        # Table VI, zhu2023: no double counting of strong/recovery slots or N_rounds/2.
        timing = dict(active_round_clock_hz=625e6, reference_fraction=160.97/625,
                      recovery_fraction=19.20/625)
        result = rate_per_second(r['bits_per_quantum_round'], timing)
        self.assertAlmostEqual(result['effective_quantum_rounds_hz'], 444.83e6, places=6)
        self.assertAlmostEqual(result['key_bps'], r['bits_per_potential_pair']*444.83e6/2, places=8)

    def test_explicit_measured_interval_sum(self):
        c = configuration()
        c['intervals'] = dict(mode='measured', delta_t_s=[1e-6, 2e-5, 1e-4], counts=[3, 7, 2])
        r = phase_average(c, .001, 128, 1.25e9)
        t = np.asarray(c['intervals']['delta_t_s'])
        # Eq.3/B10/D5, zhang2025: independent direct evaluation of the three observed bins.
        v = 4*np.pi*100*t + 4000**2*t*t
        errors = .5 - 32/(8*np.pi)*np.sin(2*np.pi/32)*np.exp(-v/2)
        self.assertAlmostEqual(r['raw_x_error'], np.dot(errors, [3, 7, 2])/12, places=14)

    def test_white_spectral_mode_matches_published_linewidth_mode(self):
        c = configuration()
        c['decoy']['error_mode'] = 'phase_weak_equal'
        c['phase']['sigma_L_rad_s'] = 0.  # Explicit numerical isolation of white-frequency contribution.
        c['intervals'] = dict(mode='measured', delta_t_s=[1e-6, 1e-5, 1e-4], counts=[2, 3, 7])
        expected = calculate(c)
        spectral = copy.deepcopy(c)
        spectral['phase']['mode'] = 'spectral'
        spectral['grid'] = dict(positive_floor_hz=1e-8, upper_hz=1e8, points=16385, interval_batch=16)
        spectral['laser_A'] = spectral['laser_B'] = dict(model='lorentzian_white_frequency', linewidth_hz=100)
        spectral['fiber_A'] = spectral['fiber_B'] = dict(length_km=0., l=0., fc1_hz=100.)
        # Eq.3, zhang2025 / Eq.1, didomenico2010: finite numerical band has a small explicit tail error.
        observed = calculate(spectral)
        self.assertLess(abs(observed['phase']['raw_x_error']-expected['phase']['raw_x_error']), 1e-6)
        self.assertLess(abs(observed['key_bps']/expected['key_bps']-1), 1e-4)

    def test_missing_phase_inputs_are_not_substituted(self):
        for key in ('sigma_L_rad_s', 'linewidth_hz', 'intrinsic_x_error'):
            c = configuration()
            del c['phase'][key]
            with self.assertRaises(KeyError):
                calculate(c)

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            pairing_rate(.1, 1.5)
        with self.assertRaises(ValueError):
            pairing_rate(1.01, 3)
        c = configuration()
        c['decoy']['QE_X'][1][1] = 2*c['decoy']['Q_X'][1][1]
        with self.assertRaisesRegex(ValueError, '<='):
            calculate(c)
        c = configuration()
        c['timing']['reference_fraction'] = 1.
        with self.assertRaisesRegex(ValueError, 'positive quantum duty'):
            calculate(c)

    def test_report_states_units_and_spool_limit(self):
        text = markdown(calculate(configuration()))
        for term in ('bit/quantum round', 'N_rounds/2', 'laboratory fiber spools',
                     'conditional residual-based 95%', 'not a single-photon error formula'):
            self.assertIn(term, text)


if __name__ == '__main__':
    unittest.main()
