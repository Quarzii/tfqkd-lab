"""Meaningful B amendment checks: diagnostic versus Eq.4, zero crossing, fit identifiability."""
import unittest
import numpy as np
from tfqkd.lab_inputs import reference_input,read_toml
from test_lab_inputs import resolve
from tfqkd.lab_engine import calculate,stability
from tfqkd.lab_diagnostics import above_band_diagnostic,balanced_arm_crossover
from tfqkd.published_fit import fit_points
from tfqkd.spectra import free_laser,free_fiber
from tfqkd.config import ROOT
import json

class AmendmentTests(unittest.TestCase):
    def test_digitized_route_matches_recorded_fit_and_exposes_limits(self):
        from tfqkd.lab_inputs import resolve as universal_resolve
        raw=read_toml('examples/bertaina2024_table3.toml')
        raw['physics'].pop('l');raw['physics'].pop('fc1_hz')
        raw['line']['spectrum']=dict(file='examples/b6/jiang2008_urban86.csv',quantity='phase',frequency_unit='Hz',psd_unit='rad^2/Hz',sidedness='one-sided',pass_='single',measurement_length_km=86)
        raw['line']['spectrum']['pass']=raw['line']['spectrum'].pop('pass_')
        lab=universal_resolve(raw)
        recorded=json.loads((ROOT/'results/stage3_b_amendments/b6/fits.json').read_text())['results'][0]['fit']
        for key,value in recorded['parameters'].items():self.assertAlmostEqual(lab.config['physics'][key]/value,1,places=6)
        self.assertTrue(any('B6_MODEL_RESIDUAL' in warning for warning in calculate(lab)['warnings']))
        self.assertFalse(lab.ranges)

    def test_above_band_does_not_mean_unreachable_threshold(self):
        lab=resolve(dict(line=dict(length_km=200,imbalance_km=0),laser=dict(preset='rio_cavity'),scheme=dict(lasers='independent',compensation='classical')))
        out=calculate(lab);d=out['classical_diagnostic']
        self.assertEqual(out['loop']['g_per_s'],0)
        self.assertGreater(d['above_band_threshold_ratio'],1)
        self.assertAlmostEqual(out['variance_rad2'],lab.config['operation']['sigma_limit_rad']**2)
        self.assertLess(out['tau_q_s'],1/d['f_b_hz'])
        self.assertEqual(d['operating_suppression_bandwidth_hz'],0)
        self.assertIn('shorter acquisition windows can',d['message'])
        manual=lab.clone();manual.settings['loop']['g_per_s']=0
        self.assertIn('manual',calculate(manual)['classical_diagnostic']['message'])

    def test_diagnostic_length_root_and_band_dependence(self):
        lab=resolve(dict(line=dict(length_km=200,imbalance_km=0),scheme=dict(lasers='independent',compensation='classical')))
        from tfqkd.config import load
        c=load();bracket=(c['classical']['length_min_km'],c['classical']['length_max_km'])
        ac=c['classical_actuator'];roots=[balanced_arm_crossover(lab,w,bracket) for w in (ac['omega_a_min_rad_s'],ac['omega_a_max_rad_s'])]
        self.assertGreater(roots[1]['arm_length_km'],roots[0]['arm_length_km'])
        for root in roots:self.assertAlmostEqual(root['ratio_at_root'],1)

    def test_fit_recovers_known_source_model(self):
        # Synthetic numerical fixture from Table III coefficients, not an equipment measurement or preset.
        from tfqkd.config import load
        c=load();p=c['physics'];f=np.geomspace(p['fc1_hz']/10,p['fc1_hz']*10,50)
        settings=json.loads((ROOT/'sources/data/b6/digitization.json').read_text())['numerics']
        fit,_=fit_points(f,free_fiber(f,114,p),dict(model='fiber',length_km=114),settings)
        for key in ('l','fc1_hz'):self.assertAlmostEqual(fit['parameters'][key]/p[key],1,places=8)
        # A flat observed frequency-noise floor cannot identify an arbitrarily distant F1 cutoff.
        f=np.geomspace(1,100,50);phase=p['r2']/f**2
        fit,_=fit_points(f,phase,dict(model='laser'),settings)
        self.assertLess(fit['jacobian_rank'],3)
        self.assertIsNone(fit['conditional_fit_se']['fc_hz'])

if __name__=='__main__':unittest.main()
