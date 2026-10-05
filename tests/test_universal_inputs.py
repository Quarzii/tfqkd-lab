"""Universal input acceptance, physics regression and inverse specification checks."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from tfqkd.config import ROOT,load
from tfqkd.lab_inputs import resolve,read_toml,MissingInputs,DEFAULT_N
from tfqkd.lab_engine import calculate,calculate_point
from tfqkd.lab_ceiling import stabilization_ceiling,table_i_ceiling_curves
from tfqkd.measured_inputs import normalize_psd,fit_psd,fit_actuator,MeasurementRejected

class UniversalTests(unittest.TestCase):
    def raw(self):return read_toml(ROOT/'examples/bertaina2024_table3.toml')

    def test_no_physical_defaults_and_complete_missing_list(self):
        with self.assertRaises(MissingInputs) as caught:resolve(dict(line=dict(length_km=200)))
        missing=caught.exception.parameters
        for path in ('line.imbalance_km','physics.fc_hz','physics.fc1_hz','keyrate.detector_efficiency','keyrate.clockrate_hz','operation.tau_ps_s','keyrate.decoy_big'):
            self.assertIn(path,missing)
        self.assertNotIn('physics.r3',missing)
        for section,key in [('laser','preset'),('line','route'),('actuator','preset'),('detector','preset')]:
            raw=self.raw();raw.setdefault(section,{})[key]='bertaina'
            with self.assertRaisesRegex(ValueError,'no equipment preset'):resolve(raw)

    def test_reference_n_is_explicitly_sourced_and_overridable(self):
        raw=self.raw();raw['physics'].pop('n');lab=resolve(raw)
        self.assertEqual(lab.config['physics']['n'],DEFAULT_N)
        row=[p for p in lab.provenance if p['parameter']=='physics.n'][0]
        self.assertTrue(row['default_used']);self.assertEqual(row['page'],2);self.assertIn('corning.com',row['source'])
        self.assertTrue(any('1.818719%' in warning for warning in lab.warnings))
        raw['physics']['n']=1.45;lab=resolve(raw)
        self.assertEqual(lab.config['physics']['n'],1.45);self.assertFalse(any('N_REFERENCE_DEFAULT' in w for w in lab.warnings))

    def test_missing_noise_generates_verified_conditional_requirement(self):
        raw=self.raw();raw['physics'].pop('r3');raw['line']['imbalance_km']=2.5
        lab=resolve(raw);self.assertEqual(lab.missing_amplitudes,('physics.r3',));self.assertEqual(lab.config['physics']['r3'],0)
        result=calculate(lab);self.assertIn('upper estimate',result['rate_estimate'])
        bound=result['equipment_requirements']['requirements'][0]
        probe=lab.clone();probe.set_parameter('physics.r3',bound['maximum'])
        actual=calculate_point(probe)
        self.assertAlmostEqual(actual['key_bps']/result['key_bps'],.9,places=6)
        self.assertGreater(bound['maximum'],0);self.assertEqual(bound['other_unmeasured_parameters_zero'],[])

    def test_zero_imbalance_does_not_claim_finite_laser_requirement(self):
        raw=self.raw();raw['physics'].pop('r3');raw['line']['imbalance_km']=0
        row=calculate(resolve(raw))['equipment_requirements']['requirements'][0]
        self.assertIsNone(row['maximum']);self.assertIn('unconstrained',row['status'])

    def test_sub_unit_and_multiple_requirements_are_not_joint_bounds(self):
        raw=self.raw()
        for name in ('r3','r2','l','s0'):raw['physics'].pop(name)
        raw['keyrate'].pop('detector_dark_count_rate_hz')
        result=calculate(resolve(raw));rows=result['equipment_requirements']['requirements']
        self.assertEqual(len(rows),5)
        for row in rows:
            self.assertAlmostEqual(row['verified_loss_fraction'],.1,places=6)
            self.assertEqual(len(row['other_unmeasured_parameters_zero']),4)
        self.assertLess(next(r for r in rows if r['parameter']=='physics.s0')['maximum'],1e-6)

    def test_passport_and_direct_inputs_agree(self):
        raw=self.raw();direct=calculate_point(resolve(raw));p=raw['physics'];k=raw['keyrate']
        raw['laser']['lorentz_width_hz']=np.pi*p.pop('r2') # Eqs.1,5 didomenico2010.
        raw['laser']['r3']=p.pop('r3');raw['laser']['fc_hz']=p.pop('fc_hz')
        raw['line']['l']=p.pop('l');raw['line']['fc1_hz']=p.pop('fc1_hz')
        raw['detector']=dict(efficiency=k.pop('detector_efficiency'),dark_count_rate_hz=k.pop('detector_dark_count_rate_hz'),error=k.pop('detector_error'))
        # User-authorized one-pole conversion: f_a=omega_a/(2*pi).
        raw['actuator']=dict(bandwidth_hz=raw['actuator']['omega_a_rad_s']/(2*np.pi))
        actual=calculate_point(resolve(raw))
        for key in ('key_bps','tau_q_s','variance_rad2'):self.assertAlmostEqual(actual[key]/direct[key],1,places=12)

    def test_psd_metadata_and_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'s.csv'
            f=np.array([1.,10.,100.]);psd=np.array([3.,4.,5.])
            np.savetxt(path,np.c_[f,psd],delimiter=',',header='frequency,psd',comments='')
            spec=dict(file=str(path),quantity='frequency',frequency_unit='Hz',psd_unit='Hz^2/Hz',sidedness='two-sided',**{'pass':'round-trip'},round_trip_psd_factor=4)
            _,converted,metadata=normalize_psd(spec)
            np.testing.assert_allclose(converted,psd*.5/f**2)
            self.assertEqual(metadata['linear_normalization_factor'],.5)
            for key in ('quantity','frequency_unit','psd_unit','sidedness','pass','round_trip_psd_factor'):
                missing=copy.deepcopy(spec);missing.pop(key)
                with self.assertRaises(MeasurementRejected):normalize_psd(missing)
            spec['psd_unit']='dBc/Hz'
            with self.assertRaises(MeasurementRejected):normalize_psd(spec)

    def test_b6_csv_acceptance_matches_original_diagnostics(self):
        settings=read_toml(ROOT/'configs/lab_defaults.toml')['fit'];historical=json.loads((ROOT/'examples/b6/historical_fits.json').read_text())['results']
        for original in historical:
            identifier=original['id'];spec=json.loads((ROOT/'examples/b6'/f'{identifier}.json').read_text())
            node='line' if original['specification']['model']=='fiber' else 'laser'
            if identifier=='jiang2008_urban86':
                parameters,metadata=fit_psd(spec,node,settings,ROOT/'examples/b6')
                for key,value in parameters.items():self.assertAlmostEqual(value/original['fit']['parameters'][key],1,places=6)
            else:
                with self.assertRaises(MeasurementRejected) as caught:fit_psd(spec,node,settings,ROOT/'examples/b6')
                metadata=caught.exception.diagnostics
            self.assertAlmostEqual(metadata['fit']['rms_log10_residual'],original['fit']['rms_log10_residual'],places=6)
            self.assertEqual(metadata['fit']['jacobian_rank'],original['fit']['jacobian_rank'])

    def test_actuator_response_identification_and_model_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'response.csv';omega=self.raw()['actuator']['omega_a_rad_s']
            # Synthetic numerical fixture from the user-adopted pole, not manufacturer data.
            f=np.geomspace(omega/(2*np.pi)/10,omega/(2*np.pi)*10,50)
            response=1/(1+2j*np.pi*f/omega)
            spec=dict(file=str(path),frequency_unit='Hz',phase_unit='deg',magnitude_unit='linear',maximum_rms_phase_deg=1.)
            settings=read_toml(ROOT/'configs/lab_defaults.toml')['fit']
            np.savetxt(path,np.c_[f,abs(response),np.angle(response,deg=True)],delimiter=',',header='frequency,magnitude,phase',comments='')
            recovered,diagnostic=fit_actuator(spec,settings)
            self.assertAlmostEqual(recovered/omega,1,places=8);self.assertIsNotNone(diagnostic['conditional_fit_se_rad_s'])
            response=response**2
            np.savetxt(path,np.c_[f,abs(response),np.angle(response,deg=True)],delimiter=',',header='frequency,magnitude,phase',comments='')
            with self.assertRaisesRegex(MeasurementRejected,'not described by one pole'):fit_actuator(spec,settings)
            np.savetxt(path,np.c_[f,np.ones_like(f),np.zeros_like(f)],delimiter=',',header='frequency,magnitude,phase',comments='')
            with self.assertRaisesRegex(MeasurementRejected,'does not identify a finite'):fit_actuator(spec,settings)

    def test_detector_psd_and_duplicate_methods_rejected(self):
        raw=self.raw();raw['detector']=dict(spectrum={})
        with self.assertRaises(ValueError):resolve(raw)
        raw=self.raw();raw['actuator']['bandwidth_hz']=1000
        with self.assertRaisesRegex(ValueError,'Choose actuator'):resolve(raw)

    def test_four_rate_ceiling_and_detector_convention(self):
        lab=resolve(self.raw());r=stabilization_ceiling(lab);rates=r['rates_bps']
        self.assertGreaterEqual(rates['R_perfect'],rates['R_classical']);self.assertGreaterEqual(rates['R_perfect'],rates['R_dual'])
        self.assertAlmostEqual(r['H'],rates['R_perfect']/rates['R_free'])
        no_det=lab.clone();no_det.set_parameter('physics.s0',0.)
        perfect=no_det.clone();perfect.settings['scheme']['compensation']='none';perfect.settings['_perfect_fiber']=True
        self.assertEqual(calculate_point(perfect)['key_bps'],rates['R_perfect'])
        self.assertTrue(r['strict_ceiling_check_passed'])

    def test_ceiling_violation_is_error(self):
        lab=resolve(self.raw())
        with patch('tfqkd.lab_ceiling.calculate_point',side_effect=[{'key_bps':1},{'key_bps':2},{'key_bps':3},{'key_bps':1}]):
            with self.assertRaisesRegex(ArithmeticError,'ceiling violated'):stabilization_ceiling(lab)

    def test_vectorized_table_i_curves_agree_with_point_core(self):
        c=load();lengths=np.array([50.,100.]);rows=table_i_ceiling_curves(c,lengths)
        from tfqkd.lab_inputs import reference_input
        for row in rows:
            sc=next(s for s in c['scenarios'] if s['name']==row['scenario'])
            for i,length in enumerate(lengths):
                lab=reference_input(c,sc,protocol=row['protocol']);lab.settings['line']=dict(arm_a_km=length+sc['delta_L_km'],arm_b_km=length)
                lab.settings['scheme']['compensation']='none'
                for name in ('free','perfect'):
                    lab.settings['_perfect_fiber']=name=='perfect'
                    point=calculate_point(lab)
                    self.assertAlmostEqual(point['key_bps']/row['R_'+name+'_bps'][i],1,places=10)

if __name__=='__main__':unittest.main()
