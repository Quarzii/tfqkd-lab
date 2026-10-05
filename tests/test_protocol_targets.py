"""Inverse protocol targets and actual-curve distance inversion, no new physics."""
import copy
import unittest
import numpy as np
from tests.gates import full_only
from tfqkd.config import ROOT
from tfqkd.lab_inputs import resolve, read_toml
from tfqkd.lab_engine import calculate, calculate_point
from tfqkd.lab_run import run
from tfqkd.lab_report import markdown, html_document
from tfqkd.lab_reach import reach_report
from tfqkd.lab_protocol_requirements import REFERENCE_F_EC


class ProtocolTargetTests(unittest.TestCase):
    def phase(self):
        return read_toml(ROOT/'examples/measured_inputs/zhou2023_403.73_ideal.toml')

    def test_protocol_inverse_recalculates_target_for_both_projections(self):
        raw = self.phase();raw['keyrate'].pop('f_error')
        # Target equals the published Table S5 rate; not an apparatus default.
        raw['requirements'] = dict(target_key_bps=146.70)
        raw['phase']['sigma_phi_rad'] = .05  # User-specified out-of-sample condition.
        lab = resolve(raw);result = run(raw, workers=1)
        text = markdown(result)
        self.assertTrue(text.startswith('Upper model estimate:'))
        self.assertIn('Upper model estimate:', html_document(text))
        rows = result['selected']['protocol_requirements']['requirements']
        self.assertEqual(len(rows), 4)
        for row in rows:
            self.assertIsNotNone(row['maximum'])
            probe = lab.clone()
            if row['parameter'] == 'keyrate.detector_error':
                probe.set_parameter('keyrate.f_error', REFERENCE_F_EC)
            probe.set_parameter(row['parameter'], row['maximum'])
            actual = calculate_point(probe)[f"projection_{row['detector_projection']}_bps"]
            self.assertLess(abs(actual/row['target_key_bps']-1), lab.config['performance']['regression_key_rtol'])
        self.assertEqual(result['reach']['status'], 'disabled')

    def test_reach_default_is_off_and_boolean_is_required(self):
        from unittest.mock import patch
        lab = resolve(read_toml(ROOT/'examples/bertaina2024_table3.toml'))
        with patch('tfqkd.lab_engine.calculate_point', side_effect=AssertionError('Disabled reach must not calculate')):
            result = reach_report(lab)
        self.assertEqual(result['status'], 'disabled')
        self.assertEqual(result['input_case_calculations'], 0)
        lab.settings['reach']['enabled'] = 'yes'
        with self.assertRaisesRegex(ValueError, 'true or false'):reach_report(lab)

    def test_missing_target_is_not_invented_and_impossible_target_is_reported(self):
        raw = self.phase();result = calculate(resolve(raw))
        self.assertIsNone(result['protocol_requirements']['target_key_bps'])
        raw['requirements'] = dict(target_key_bps=2*result['key_bps'])  # Numerical unattainable-target test.
        rows = calculate(resolve(raw))['protocol_requirements']['requirements']
        self.assertTrue(all(r['maximum'] is None and 'unattainable' in r['status'] for r in rows))
        raw['requirements']['target_key_bps'] = np.nan
        with self.assertRaises(ValueError):resolve(raw)

    def test_reach_prunes_only_positive_monotone_envelopes(self):
        from unittest.mock import patch
        from tfqkd.lab_reach import _length_point
        lab=resolve(read_toml(ROOT/'examples/bertaina2024_table3.toml'))
        coordinates=[{'physics.l':44.},{'physics.l':88.}] # Table III base and explicit numerical factor-two stress.
        # Numerical scheduling fixtures; signed negative rates are deliberately not monotone.
        positive=dict(key_bps=1.,raw_key_bps=1.,spectral_calculations=1)
        with patch('tfqkd.lab_engine.calculate_point',return_value=positive) as calculate:
            point,cases,_=_length_point(lab,coordinates,.02,200.,('physics.l',))
            self.assertEqual(cases,1);self.assertEqual(calculate.call_args.args[0].config['physics']['l'],88.)
        negative=[dict(key_bps=0.,raw_key_bps=-1.,spectral_calculations=1),
                  dict(key_bps=0.,raw_key_bps=-2.,spectral_calculations=1)]
        with patch('tfqkd.lab_engine.calculate_point',side_effect=negative):
            point,cases,_=_length_point(lab,coordinates,.02,200.,('physics.l',))
            self.assertEqual(cases,2);self.assertEqual(point['raw_key_bps'],-2.)

    def test_spectral_inputs_allow_only_labelled_ideal_protocol_values(self):
        raw = read_toml(ROOT/'examples/bertaina2024_table3.toml')
        raw['keyrate'].pop('detector_error');raw['keyrate'].pop('f_error')
        lab = resolve(raw)
        self.assertEqual(set(lab.settings['_ideal_protocol_fields']), {'detector_error', 'f_error'})
        result = calculate(lab)
        self.assertIn('upper model estimate', result['rate_estimate'])
        lab.set_parameter('keyrate.detector_error', .02)  # Table II published reference fixture.
        self.assertNotIn('detector_error', lab.settings['_ideal_protocol_fields'])

    @full_only
    def test_reach_inverts_actual_curve_and_verifies_factor_two(self):
        raw = read_toml(ROOT/'examples/bertaina2024_table3.toml')
        probe = copy.deepcopy(raw);probe['line']['length_km'] = 400.
        target = calculate_point(resolve(probe))['key_bps']
        self.assertGreater(target, 0)
        # Explicit numerical validation domain; never an installed-line prior.
        raw['reach'] = dict(enabled=True, minimum_total_km=100., maximum_total_km=1000., geometry='fixed_imbalance', points=17)
        raw['requirements'] = dict(target_key_bps=target)
        lab = resolve(raw);reach = reach_report(lab)
        serial = reach_report(lab, workers=1)
        self.assertEqual(reach['zero_rate_reach_km'], serial['zero_rate_reach_km'])
        self.assertEqual(reach['target_reach'], serial['target_reach'])
        self.assertAlmostEqual(reach['target_reach']['nominal']['reach_km'], 400., places=3)
        for name, multiplier in [('nominal', 1.), ('doubled_rate', .5), ('halved_rate', 2.)]:
            length = reach['target_reach'][name]['reach_km']
            candidate = lab.clone();candidate.settings['line']['length_km'] = length
            actual = calculate_point(candidate)['key_bps']
            self.assertLess(abs(actual/(target*multiplier)-1), lab.config['performance']['regression_key_rtol'])
        self.assertGreater(reach['doubling_shift_km'], 0)
        self.assertGreater(reach['halving_shift_km'], 0)
        self.assertIsNotNone(reach['zero_rate_reach_km'])
        self.assertIn('A factor of 2 in key rate corresponds to', reach['message'])
        automatic = copy.deepcopy(raw);automatic['reach'] = dict(enabled=True)
        auto = reach_report(resolve(automatic))
        self.assertIn('automatic', auto['domain_source'])
        self.assertAlmostEqual(auto['zero_rate_reach_km'], reach['zero_rate_reach_km'], places=3)
        self.assertAlmostEqual(auto['doubling_shift_km'], reach['doubling_shift_km'], places=3)
        raw['reach']['maximum_total_km'] = 300.
        censored = reach_report(resolve(raw))
        self.assertIsNone(censored['target_reach']['nominal']['reach_km'])
        self.assertIn('right-censored', censored['target_reach']['nominal']['status'])

    @full_only
    def test_reach_uses_lower_existing_range_envelope_at_each_length(self):
        raw = read_toml(ROOT/'examples/bertaina2024_table3.toml')
        nominal = raw['physics']['l']
        # Numerical factor-two stress around Table III; NOT an empirical equipment range.
        raw['ranges'] = {'physics.l':dict(minimum=nominal, maximum=2*nominal, sources=['Table III Bertaina coefficient; explicitly labelled numerical factor-two stress'])}
        raw['reach'] = dict(enabled=True, minimum_total_km=100., maximum_total_km=1000., geometry='fixed_imbalance', points=17)
        probe = copy.deepcopy(raw);probe.pop('ranges');probe['line']['length_km'] = 400.
        target = calculate_point(resolve(probe))['key_bps']
        raw['requirements'] = dict(target_key_bps=target)
        cases = [dict(coordinates={'physics.l':value}) for value in (nominal, 2*nominal)]
        reach = reach_report(resolve(raw), range_cases=cases)
        scheduled=reach_report(resolve(raw),range_cases=cases,monotone_parameters=('physics.l',))
        self.assertEqual(reach['target_reach'],scheduled['target_reach'])
        self.assertEqual(reach['zero_rate_reach_km'],scheduled['zero_rate_reach_km'])
        self.assertIn('lower sampled', reach['rate_curve_basis'])
        length = reach['target_reach']['nominal']['reach_km']
        rates = []
        for value in (nominal, 2*nominal):
            point = copy.deepcopy(probe);point['physics']['l'] = value;point['line']['length_km'] = length
            rates.append(calculate_point(resolve(point))['key_bps'])
        self.assertLess(abs(min(rates)/target-1), resolve(raw).config['performance']['regression_key_rtol'])


if __name__ == '__main__':unittest.main()
