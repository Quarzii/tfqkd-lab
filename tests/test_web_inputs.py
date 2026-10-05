"""Web contract checks: range errors, normalization and published examples."""
import copy
import json
import math
import tomllib
import unittest
from unittest import mock
from pathlib import Path
from tfqkd.config import ROOT
from tfqkd.web_inputs import FieldValidationError, validate_web_input
from tfqkd.lab_web import _example, _canonical_toml, _web_context
from tfqkd.lab_inputs import resolve
from tfqkd.lab_engine import calculate_point
from tests.test_lab_web import request

class WebInputTests(unittest.TestCase):
    def setUp(self):
        self.raw=tomllib.loads((ROOT/'examples/measured_inputs/zhou2023_403.73_ideal.toml').read_text())

    def test_boundaries_with_specific_field_errors_before_calculation(self):
        cases=[('phase','sigma_phi_rad',math.pi+1e-6),('phase','sigma_phi_rad',-.01),
               ('keyrate','pz_sns',0),('keyrate','eps_sns_aopp',1.01),
               ('keyrate','f_error',.99),('keyrate','clockrate_hz',0),
               ('keyrate','detector_error',.5),('line','loss_a_db',0),('line','arm_a_km',-1)]
        for section,key,value in cases:
            with self.subTest(section=section,key=key):
                raw=copy.deepcopy(self.raw);raw[section][key]=value
                with mock.patch('tfqkd.lab_run.run') as run:
                    status,_,body=request('POST','/api/run',{'toml':_canonical_toml(raw)})
                self.assertEqual(status,422)
                self.assertIn(section+'.'+key,json.loads(body)['field_errors']);run.assert_not_called()
        for efficiency in (0,1.01):
            raw=copy.deepcopy(self.raw);raw['detector']['channels'][0]['efficiency']=efficiency
            with self.assertRaises(FieldValidationError) as caught:validate_web_input(raw)
            self.assertIn('detector.channels.0.efficiency',caught.exception.field_errors)
        raw=copy.deepcopy(self.raw);raw['detector']['channels'][1]['dark_count_rate_hz']=-1
        with self.assertRaises(FieldValidationError):validate_web_input(raw)

    def test_allowed_boundaries_and_unusual_values(self):
        raw=copy.deepcopy(self.raw)
        raw['phase']['sigma_phi_rad']=math.pi
        raw['detector']['channels'][0]['efficiency']=1
        raw['detector']['channels'][0]['dark_count_rate_hz']=0
        raw['keyrate'].update(f_error=1,detector_error=0,eps_sns_aopp=.269)
        validate_web_input(raw)
        raw['phase']['sigma_phi_rad']=1.2
        raw['detector']['channels'][0]['efficiency']=.04
        raw['keyrate']['f_error']=2.1
        warnings=validate_web_input(raw)
        self.assertEqual(len(warnings),3)
        status,_,body=request('POST','/api/run',{'toml':_canonical_toml(raw)})
        self.assertEqual(status,200)
        self.assertEqual(len(json.loads(body)['result']['web_input']['warnings']),3)

    def test_epsilon_requires_both_event_types(self):
        for epsilon in (0,1):
            raw=copy.deepcopy(self.raw);raw['keyrate']['eps_sns_aopp']=epsilon
            with self.assertRaises(FieldValidationError) as caught:validate_web_input(raw)
            self.assertIn('both sending and not-sending',caught.exception.field_errors['keyrate.eps_sns_aopp'])
        for epsilon in (.01,.5):
            raw=copy.deepcopy(self.raw);raw['keyrate']['eps_sns_aopp']=epsilon
            self.assertIn('keyrate.eps_sns_aopp',validate_web_input(raw))

    def test_decoy_order_reports_all_three_fields(self):
        self.raw['keyrate']['decoy_big']=self.raw['keyrate']['decoy_medium']
        with self.assertRaises(FieldValidationError) as caught:validate_web_input(self.raw)
        self.assertEqual(set(caught.exception.field_errors),{'keyrate.decoy_big','keyrate.decoy_medium','keyrate.decoy_mini'})

    def test_losses_override_attenuation_even_in_advanced(self):
        self.raw['line']['attenuation_db_per_km']=0
        self.raw['keyrate']['attenuation_db_per_km']=.2
        original=_canonical_toml(self.raw)
        context=_web_context(self.raw,{})
        self.assertNotIn('attenuation_db_per_km',self.raw['line'])
        self.assertNotIn('attenuation_db_per_km',self.raw['keyrate'])
        self.assertEqual(len(context['ignored_attenuation']),2)
        status,_,body=request('POST','/api/run',{'toml':original})
        self.assertEqual(status,200)
        self.assertIn('per-km attenuation was not used',json.loads(body)['markdown'])
        self.assertEqual(len(json.loads(body)['result']['web_input']['ignored_attenuation']),2)

    def test_zhou_per_user_input_preserves_held_out_predictions(self):
        evidence=json.loads((ROOT/'results/out_of_sample/validation.json').read_text())
        per_user={'decoy_big':.493,'decoy_medium':.105,'decoy_mini':.0002}
        for row in evidence['zhou']:
            anchor=str(row['anchor_km'])
            raw=tomllib.loads((ROOT/f'examples/measured_inputs/zhou2023_{anchor}_ideal.toml').read_text())
            raw['phase']['sigma_phi_rad']=row['sigma_rad'];raw['keyrate']['detector_error']=row['e_d']
            # Appendix D, bertaina2024 / author Cell21: GUI normalization to total SNS intensities.
            # Different lengths have different source intensities; use each article point's values.
            published=tomllib.loads((ROOT/f'examples/measured_inputs/zhou2023_{anchor}_ideal.toml').read_text())
            point_per_user={k:published['keyrate'][k]/2 for k in per_user}
            raw['keyrate'].update({k:2*v for k,v in point_per_user.items()})
            context=_web_context(raw,{'intensities_per_user':point_per_user})
            own=calculate_point(resolve(raw));projection=int(row['detector'][1:])
            self.assertAlmostEqual(own[f'projection_{projection}_bps'],row['fitted_rate_bps'],places=9)
            self.assertEqual([r['factor'] for r in context['intensities']],[2,2,2])
            for prediction in row['predictions']:
                length=str(prediction['length_km'])
                target=tomllib.loads((ROOT/f'examples/measured_inputs/zhou2023_{length}_ideal.toml').read_text())
                displayed={k:target['keyrate'][k]/2 for k in per_user}
                target['keyrate'].update({k:2*v for k,v in displayed.items()})
                target['phase']['sigma_phi_rad']=row['sigma_rad'];target['keyrate']['detector_error']=row['e_d']
                _web_context(target,{'intensities_per_user':displayed})
                result=calculate_point(resolve(target))
                self.assertAlmostEqual(result[f'projection_{projection}_bps'],prediction['predicted_key_bps'],places=9)

    def test_examples_and_report_show_both_intensities(self):
        for name in ('zhou2023','bertaina2024'):
            status,_,body=request('GET','/api/examples/'+name)
            self.assertEqual(status,200);self.assertIn('source',json.loads(body))
        b=_example('bertaina2024')['fields']
        self.assertEqual((b['laserModel'],b['schemeLasers'],b['schemeCompensation']),('free','common','dual'))
        loaded=_example('zhou2023')
        self.assertEqual([loaded['fields'][key] for key in ('decoyBig','decoyMedium','decoyMini')],[.493,.105,.0002])
        status,_,body=request('POST','/api/run',{'toml':loaded['configuration'],'intensities_per_user':{'decoy_big':.493,'decoy_medium':.105,'decoy_mini':.0002},'example_source':loaded['source']})
        self.assertEqual(status,200);payload=json.loads(body)
        self.assertIn('Per user [photons/pulse]',payload['markdown'])
        self.assertIn('0.493',payload['markdown']);self.assertIn('0.986',payload['markdown'])
        with self.assertRaises(FieldValidationError):_web_context(copy.deepcopy(self.raw),{'intensities_per_user':{'decoy_big':.986}})

    def test_cal_intensity_is_not_doubled(self):
        raw={'protocol':{'name':'CAL'},'keyrate':{'u_cal':.018}}
        context=_web_context(raw,{'intensities_per_user':{'u_cal':.018}})
        self.assertEqual(context['intensities'][0]['core'],.018)
        self.assertEqual(context['intensities'][0]['factor'],1)

if __name__=='__main__':unittest.main()
