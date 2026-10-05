"""Independent checks of the authorized actuator pole and adaptive key-rate criterion."""
import unittest
import numpy as np
from scipy.optimize import brentq
from tfqkd.config import load
from tfqkd.classical import (delay_and_boundary, actuator_stability, actuator_loop_gain,
                            actuator_controller, actuator_pole, require_stable_actuator,
                            remote_power_ratio, remote_controller_ratio)
from tfqkd.classical_keyrate import converged_operating_point
from tfqkd.keyrates import KeyRateParameters, sns_aopp_rate_per_pulse, cal_rate_per_pulse


class ActuatorTests(unittest.TestCase):
    def test_critical_gain_and_root_direction(self):
        c = load(); ac = c['classical_actuator']
        _, rt, _ = delay_and_boundary(c['classical']['length_km'], c['physics'])
        for wa in [ac['omega_a_min_rad_s'], ac['omega_a_max_rad_s']]:
            st = actuator_stability(rt, wa)
            # Eq. A6 with the adopted pole: verify the original complex loop, not its reduced real equation.
            g_loop = actuator_loop_gain([st['critical_frequency_hz']], rt, st['g_crit_per_s'], wa)
            np.testing.assert_allclose(g_loop, -1, atol=c['validation']['reference_rtol'])
            for fraction in [1 - ac['root_crossing_fraction'], 1 + ac['root_crossing_fraction']]:
                z = actuator_pole(fraction * st['g_crit_per_s'] * rt, st['b'], 1j * st['x'])
                self.assertEqual(z.real > 0, fraction > 1)

    def test_first_crossing_is_smallest_gain(self):
        c = load(); ac = c['classical_actuator']
        _, rt, _ = delay_and_boundary(c['classical']['length_km'], c['physics'])
        for wa in [ac['omega_a_min_rad_s'], ac['omega_a_max_rad_s']]:
            gains = [actuator_stability(rt, wa, index)['g_crit_per_s'] for index in range(c['classical']['pole_branches'])]
            self.assertTrue(np.all(np.diff(gains) > 0))

    def test_unstable_and_marginal_gain_rejected(self):
        c = load(); ac = c['classical_actuator']
        _, rt, _ = delay_and_boundary(c['classical']['length_km'], c['physics'])
        wa = ac['omega_a_min_rad_s']; gc = actuator_stability(rt, wa)['g_crit_per_s']
        for g in [gc, gc * (1 + ac['root_crossing_fraction'])]:
            with self.assertRaises(ValueError):
                require_stable_actuator(rt, g, wa)

    def test_remote_integral_with_actuator_against_spatial_quadrature(self):
        c = load(); ac = c['classical_actuator']
        tau, rt, _ = delay_and_boundary(c['classical']['length_km'], c['physics'])
        for wa in [ac['omega_a_min_rad_s'], ac['omega_a_max_rad_s']]:
            st = actuator_stability(rt, wa)
            controller = actuator_controller(ac['gain_fractions'][0] * st['g_crit_per_s'], wa)
            f = np.geomspace(1 / c['operation']['tau_max_s'], st['critical_frequency_hz'] * ac['quadrature_frequency_multiple'], ac['quadrature_frequency_points'])
            expected = remote_power_ratio(f, tau, controller, c['classical']['spatial_nodes'])
            actual = remote_controller_ratio(f, tau, controller)
            np.testing.assert_allclose(actual, expected, rtol=c['validation']['reference_rtol'])

    def test_operating_time_is_solved_and_not_fixed_at_cap(self):
        c = load(); ac = c['classical_actuator']; length = c['classical']['length_km']
        _, rt, _ = delay_and_boundary(length, c['physics'])
        wa = ac['omega_a_min_rad_s']; g = ac['gain_fractions'][0] * actuator_stability(rt, wa)['g_crit_per_s']
        row, _ = converged_operating_point(c, length, 'classical', g, wa)
        self.assertEqual(row['operating_status'], 'threshold')
        self.assertLess(row['tau_q_s'], c['operation']['tau_max_s'])
        self.assertAlmostEqual(row['sigma_phi_rad'], c['operation']['sigma_limit_rad'])
        dual, _ = converged_operating_point(c, length, 'dual')
        self.assertEqual(dual['operating_status'], 'capped')
        self.assertEqual(dual['tau_q_s'], c['operation']['tau_max_s'])

    def test_key_zero_is_not_fixed_time_phase_crossing(self):
        c = load(); ac = c['classical_actuator']; params = KeyRateParameters(**c['keyrate'])
        for protocol, model in [('sns', sns_aopp_rate_per_pulse), ('cal', cal_rate_per_pulse)]:
            # Appendix A Eq. A1 / Appendix B Eq. B1: independent key-bound zero at sigma threshold.
            loss = brentq(lambda value: model(value, c['operation']['sigma_limit_rad'], params),
                          ac['key_zero_loss_min_db'], ac['key_zero_loss_max_db'])
            length = loss / (2 * params.attenuation_db_per_km)
            _, rt, _ = delay_and_boundary(length, c['physics'])
            wa = ac['omega_a_min_rad_s']; g = ac['gain_fractions'][0] * actuator_stability(rt, wa)['g_crit_per_s']
            row, _ = converged_operating_point(c, length, 'classical', g, wa)
            self.assertEqual(row['operating_status'], 'threshold')
            self.assertLess(row['tau_q_s'], c['operation']['tau_max_s'])
            self.assertAlmostEqual(row[protocol + '_raw_bps'], 0, places=7)


if __name__ == '__main__':
    unittest.main()
