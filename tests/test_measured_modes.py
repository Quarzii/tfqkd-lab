"""Meaningful checks of measured wrappers against unchanged spectral/protocol kernels."""
import copy
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from tests.gates import full_only
from unittest.mock import patch
import numpy as np
from tfqkd.config import ROOT,load
from tfqkd.lab_inputs import resolve,read_toml,MissingInputs,reference_input
from tfqkd.lab_engine import calculate_point
from tfqkd.lab_run import run
from tfqkd.lab_report import write_report
from tfqkd.keyrates import KeyRateParameters,sns_aopp_rate_per_pulse,cal_rate_per_pulse
from tfqkd.integration import frequency_grid
from tfqkd.spectra import free_laser,stabilized_laser,free_fiber,stabilized_fiber
from tfqkd.measured_inputs import MeasurementRejected


class MeasuredModesTests(unittest.TestCase):
    def raw(self):return read_toml(ROOT/'examples/bertaina2024_table3.toml')

    def phase_input(self):
        raw=self.raw();point=calculate_point(resolve(raw))
        # Executed reference fixture, not invented measurement coefficients.
        return dict(phase=dict(mode='measured',sigma_phi_rad=point['sigma_phi_rad'],tau_s=point['tau_q_s'],tau_ps_s=raw['operation']['tau_ps_s']),
                    line=copy.deepcopy(raw['line']),keyrate=copy.deepcopy(raw['keyrate']),protocol=raw['protocol'])

    def test_phase_bypasses_spectral_chain_and_matches_kernel(self):
        raw=self.raw();expected=calculate_point(resolve(raw));observed=self.phase_input()
        with patch('tfqkd.lab_engine.PhaseIntegral',side_effect=AssertionError('spectral chain used')),patch('tfqkd.lab_engine.frequency_grid',side_effect=AssertionError('frequency grid used')):
            result=run(observed,workers=1)
        for key in ('key_bps','variance_rad2','tau_q_s','duty'):self.assertAlmostEqual(result['selected'][key]/expected[key],1,places=12)
        self.assertEqual(result['spectral_calculations'],0)
        self.assertNotIn('fiber_ceiling',result);self.assertNotIn('variance_contributions',result)
        np.testing.assert_allclose(result['selected']['e_phi'],expected['variance_rad2']/4,rtol=load()['validation']['reference_rtol']) # Eq.1, bertaina2024; sqrt/square roundoff is permitted.
        with tempfile.TemporaryDirectory() as directory:
            write_report(result,directory)
            self.assertFalse(list(Path(directory).glob('*.png')))
            text=(Path(directory)/'report.md').read_text();html=(Path(directory)/'report.html').read_text()
            self.assertIn('not a spectral forecast',text);self.assertIn('no phase-threshold solve',html)

    def test_phase_required_inputs_and_ideal_errors_are_explicit(self):
        raw=self.phase_input()
        for field in ('sigma_phi_rad','tau_s','tau_ps_s'):
            probe=copy.deepcopy(raw);probe['phase'].pop(field)
            with self.assertRaises(MissingInputs):resolve(probe)
        raw['keyrate'].pop('detector_error');raw['keyrate'].pop('f_error')
        result=run(raw,workers=1)
        self.assertIn('upper model estimate',result['selected']['rate_estimate'])
        self.assertTrue(any('IDEAL_PROTOCOL' in w for w in result['selected']['warnings']))

    def test_total_losses_and_two_real_detector_projections(self):
        raw=self.phase_input();params=KeyRateParameters(**raw['keyrate'])
        a,b=resolve(self.raw()).arms()
        # QKD.ipynb loss convention; explicit losses from the same Table II fixture.
        raw['line'].update(loss_a_db=a*params.attenuation_db_per_km,loss_b_db=b*params.attenuation_db_per_km)
        raw['keyrate'].pop('attenuation_db_per_km');raw['keyrate'].pop('detector_efficiency');raw['keyrate'].pop('detector_dark_count_rate_hz')
        raw['detector']=dict(channels=[dict(name='D0',efficiency=.60,dark_count_rate_hz=2.,background_count_rate_hz=2.),dict(name='D1',efficiency=.65,dark_count_rate_hz=2.,background_count_rate_hz=2.)]) # Zhou2023 Note3/Table S3; numerical receiver fixture, not this Bertaina apparatus.
        result=run(raw,workers=1)['selected'];loss=2*max(raw['line']['loss_a_db'],raw['line']['loss_b_db'])
        for channel in result['receiver']['channels']:
            p=replace(params,detector_efficiency=channel['efficiency'],detector_dark_count_rate_hz=channel['protocol_noise_count_rate_hz'])
            # Eq.A1, bertaina2024; unchanged scalar protocol.
            expected=p.clockrate_hz*result['duty']*sns_aopp_rate_per_pulse(loss,result['sigma_phi_rad'],p)
            self.assertAlmostEqual(channel['raw_key_bps']/expected,1,places=12)
        self.assertEqual(result['key_bps'],min(r['key_bps'] for r in result['receiver']['channels']))
        self.assertIn('No detector averaging',result['receiver']['reduction'])
        # Without geometry, complete measured losses still suffice in the phase-observation mode.
        raw['line']={k:v for k,v in raw['line'].items() if k.startswith('loss_')}
        self.assertIsNone(run(raw,workers=1)['selected']['arm_lengths_km'])

    def spectrum(self,path):
        return dict(file=str(path),mode='direct',quantity='phase',frequency_unit='Hz',psd_unit='rad^2/Hz',sidedness='one-sided',**{'pass':'single'})

    @full_only
    def test_direct_all_table_i_output_spectra_match_parameter_path(self):
        c=load();c['grid']['mode']='fast';f=frequency_grid(c,mode='reference')
        with tempfile.TemporaryDirectory() as directory:
            for sc in c['scenarios']:
                lab=reference_input(c,sc);expected=calculate_point(lab);raw=self.raw()
                raw['physics']=copy.deepcopy(lab.config['physics']);raw['line']=copy.deepcopy(lab.settings['line']);raw['scheme']=copy.deepcopy(lab.settings['scheme']);raw['laser']=copy.deepcopy(lab.settings['laser'])
                laser=(stabilized_laser if sc['cavity'] else free_laser)(f,c['physics']) # F1–F4, bertaina2024.
                line=(stabilized_fiber if sc['stabilized'] else free_fiber)(f,lab.arms()[1],c['physics']) # Eq.6/8, bertaina2024; detection excluded.
                for node,psd in [('laser',laser),('line',line)]:
                    path=Path(directory)/(node+'.csv');np.savetxt(path,np.c_[f,psd],delimiter=',',header='frequency,psd',comments='')
                    raw[node]['spectrum']=self.spectrum(path)
                raw['line']['spectrum']['measurement_length_km']=lab.arms()[1]
                for key in ('r3','r2','fc_hz','C4','C3','C2','B_hz','gamma','delta','l','fc1_hz','lambda_s_nm','lambda_q_nm'):raw['physics'].pop(key,None)
                actual=calculate_point(resolve(raw))
                for key in ('variance_rad2','tau_q_s','duty','key_bps'):
                    self.assertLess(abs(actual[key]/expected[key]-1),c['performance']['regression_variance_rtol'])

    @full_only
    def test_direct_band_rejection_extrapolation_and_analysis_omission(self):
        raw=self.raw();raw['scheme']['compensation']='none';raw['line'].pop('imbalance_km');raw['line']['arm_a_km']=raw['line'].pop('length_km')/2;raw['line']['arm_b_km']=raw['line']['arm_a_km']
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'laser.csv';f=np.geomspace(1.,1e3,50) # Numerical truncated-band fixture.
            np.savetxt(path,np.c_[f,free_laser(f,raw['physics'])],delimiter=',',header='frequency,psd',comments='')
            raw['laser']['spectrum']=self.spectrum(path)
            for key in ('r3','r2','fc_hz','C4','C3','C2','B_hz','gamma','delta'):raw['physics'].pop(key,None)
            with self.assertRaisesRegex(MeasurementRejected,'exceeds measured band'):resolve(raw)
            raw['laser']['spectrum']['extrapolation']='power-law'
            result=run(raw,workers=1)
            self.assertTrue(any('EXTRAPOLATION' in w for w in result['selected']['warnings']))
            self.assertTrue(any(r['parameter']=='physics.r3' and 'direct PSD' in r['status'] for r in result['sensitivity']['not_ranked']))
            self.assertFalse(any(r['parameter'].startswith('physics.r') for r in result['sensitivity']['ranked']))
            self.assertNotIn('physics.r3',result['selected']['missing_amplitudes'])
            # End-to-end export must serialize normalized CSV metadata, including boolean extrapolation flags.
            report_dir=Path(directory)/'report';write_report(result,report_dir)
            self.assertIn('EXTRAPOLATION',(report_dir/'report.md').read_text())
            self.assertTrue((report_dir/'result.json').is_file())

    def test_identical_separate_detectors_and_losses_preserve_scalar_path(self):
        raw=self.raw();lab=resolve(raw);expected=calculate_point(lab);a,b=lab.arms();eta=raw['keyrate'].pop('detector_efficiency');dark=raw['keyrate'].pop('detector_dark_count_rate_hz');alpha=raw['keyrate'].pop('attenuation_db_per_km')
        raw['line'].update(loss_a_db=a*alpha,loss_b_db=b*alpha)
        raw['detector']=dict(channels=[dict(name=name,efficiency=eta,dark_count_rate_hz=dark) for name in ('D0','D1')])
        actual=calculate_point(resolve(raw))
        for key in ('key_bps','raw_key_bps','tau_q_s','variance_rad2','duty'):self.assertEqual(actual[key],expected[key])


if __name__=='__main__':unittest.main()
