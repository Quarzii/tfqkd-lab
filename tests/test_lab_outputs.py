"""Part C checks: applied upgrades, omissions, additivity and report boundaries."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from tests.gates import full_only
import numpy as np
from tfqkd.config import ROOT,load
from tfqkd.lab_inputs import resolve,read_toml,reference_input
from tfqkd.lab_engine import calculate_point
from tfqkd.lab_analysis import sensitivity,variance_contributions
from tfqkd.lab_limits import tool_limitations_markdown,limitation_ids
from tfqkd.lab_run import run,markdown
from tfqkd.lab_report import html_document,write_report
from tfqkd.measured_inputs import fit_actuator,MeasurementRejected
from tfqkd.lab_compare import compare,comparison_markdown

class OutputTests(unittest.TestCase):
    def raw(self):return read_toml(ROOT/'examples/bertaina2024_table3.toml')

    def test_shortfall_transform_and_clock_explanation(self):
        lab=resolve(self.raw());point=calculate_point(lab);rows=sensitivity(lab,point,workers=1)['ranked']
        eff=next(r for r in rows if r['parameter']=='keyrate.detector_efficiency')
        self.assertEqual(eff['changed_value'],.95);self.assertEqual(eff['operation'],'halve shortfall to unity')
        clock=next(r for r in rows if r['parameter']=='keyrate.clockrate_hz')
        self.assertAlmostEqual(clock['key_bps']/point['key_bps'],2,delta=.01)
        self.assertIn('not due to noise suppression',clock['note'])
        self.assertEqual([r['gain_bps'] for r in rows],sorted([r['gain_bps'] for r in rows],reverse=True))
        for excluded in ('physics.n','physics.fc_hz','physics.fc1_hz','physics.fc2_hz','physics.lambda_s_nm','keyrate.eps_sns_aopp'):
            self.assertNotIn(excluded,[r['parameter'] for r in rows])

    def test_unknown_amplitude_not_used_as_sensitivity(self):
        raw=self.raw();raw['physics'].pop('r3');raw['line']['imbalance_km']=2.5
        lab=resolve(raw);point=calculate_point(lab);result=sensitivity(lab,point,workers=1)
        self.assertNotIn('physics.r3',[r['parameter'] for r in result['ranked']])
        missing=next(r for r in result['not_ranked'] if r['parameter']=='physics.r3')
        self.assertIn('equipment requirement',missing['status'])

    def test_variance_additivity_for_table_i(self):
        c=load()
        for sc in c['scenarios']:
            lab=reference_input(c,sc);point=calculate_point(lab);d=variance_contributions(lab,point)
            self.assertLess(d['sum_relative_error'],1e-12)
            self.assertAlmostEqual(sum(r['fraction'] for r in d['rows']),1)
            detection=next(r for r in d['rows'] if r['component']=='detection')['variance_rad2']
            self.assertEqual(detection>0,sc['stabilized'])

    def test_zero_variance_has_no_invented_fraction(self):
        raw=self.raw()
        for key in ('r3','r2','l','s0'):raw['physics'][key]=0.
        lab=resolve(raw);point=calculate_point(lab);d=variance_contributions(lab,point)
        self.assertEqual(d['total_variance_rad2'],0)
        self.assertIsNone(d['dominant_component'])
        self.assertTrue(all(row['fraction'] is None for row in d['rows']))

    def test_classical_variance_uses_winner_adaptive_grid(self):
        raw=self.raw();raw['laser']['model']='cavity';raw['scheme']=dict(lasers='independent',compensation='classical')
        raw['line']=dict(length_km=100.,imbalance_km=0.)
        raw['loop']=dict(g_per_s=100.) # Numerical stable-g fixture; not an apparatus default.
        lab=resolve(raw);point=calculate_point(lab);d=variance_contributions(lab,point)
        self.assertEqual(len(d['frequency_hz']),point['grid_points']);self.assertLess(d['sum_relative_error'],1e-12)
        self.assertEqual(next(r for r in d['rows'] if r['component']=='detection')['variance_rad2'],0)
        self.assertTrue(any('CLASSICAL_DETECTION_NOISE' in w for w in point['warnings']))

    def test_explicit_arm_imbalance_sensitivity_keeps_total(self):
        raw=self.raw();raw['line']=dict(arm_a_km=101.25,arm_b_km=98.75)
        lab=resolve(raw);point=calculate_point(lab);rows=sensitivity(lab,point,workers=1)['ranked']
        imbalance=next(r for r in rows if r['parameter']=='line.imbalance_km')
        self.assertEqual(imbalance['changed_value'],1.25)
        probe=lab.clone();probe.settings['line']=dict(length_km=200.,imbalance_km=1.25)
        self.assertEqual(imbalance['key_bps'],calculate_point(probe)['key_bps'])

    def test_phase_fit_default_five_degrees_and_override(self):
        settings=read_toml(ROOT/'configs/lab_defaults.toml')['fit'];omega=self.raw()['actuator']['omega_a_rad_s']
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'response.csv';f=np.geomspace(100,10000,50)
            # Synthetic user-authorized one-pole fit fixture with a 3-degree phase offset; not measured hardware data.
            response=1/(1+2j*np.pi*f/omega)
            np.savetxt(path,np.c_[f,abs(response),np.angle(response,deg=True)+3.],delimiter=',',header='frequency,magnitude,phase',comments='')
            spec=dict(file=str(path),frequency_unit='Hz',phase_unit='deg',magnitude_unit='linear')
            fitted,metadata=fit_actuator(spec,settings)
            self.assertEqual(metadata['acceptance_phase_limit_deg'],5)
            self.assertTrue(metadata['acceptance_phase_default_used'])
            spec['maximum_rms_phase_deg']=1.
            with self.assertRaises(MeasurementRejected):fit_actuator(spec,settings)
            spec.pop('maximum_rms_phase_deg');settings=dict(settings,maximum_rms_phase_deg=1.)
            with self.assertRaises(MeasurementRejected):fit_actuator(spec,settings)

    def test_report_order_and_only_tool_limitations(self):
        result=run(self.raw(),workers=1);text=markdown(result);document=html_document(text)
        headers=['1. Fiber stabilization ceiling','2. Selected installation','3. Ranked sensitivity','4. Conditional equipment requirements','5. Recommended classical loop','6. Phase-variance contributions','7. Applicability and Tool limitations']
        self.assertEqual([text.index(h) for h in headers],sorted(text.index(h) for h in headers))
        self.assertIn('CLASSICAL_DETECTION_NOISE',text);self.assertIn('optimistic',text)
        for identifier in ('T4_UNCERTAINTY','T4_HIGH_FREQUENCY','B6_RMS_CLOSURE','APPARATUS_IDENTITY','Research-stage assumptions','Validation limits'):
            self.assertNotIn(identifier,text);self.assertNotIn(identifier,document)
        self.assertIn('CONDITIONAL_NOISE_REQUIREMENTS',result['limits'])
        self.assertNotIn('T4_UNCERTAINTY',result['limits'])
        malicious=html_document('# <script>alert(1)</script>\n\n| Name | Value |\n| --- | --- |\n| <img src=x onerror=alert(1)> | 1 |')
        self.assertNotIn('<script>',malicious);self.assertNotIn('<img src=x',malicious)

    @full_only
    def test_two_and_three_variants_and_working_ratios(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory);example=(ROOT/'examples/bertaina2024_table3.toml').read_text()
            (base/'dual.toml').write_text(example)
            (base/'none.toml').write_text(example.replace('compensation = "dual"','compensation = "none"'))
            (base/'classical.toml').write_text(example.replace('compensation = "dual"','compensation = "classical"')+'\n[loop]\ng_per_s = 0.0 # Explicit test, not instrument default\n')
            spec=dict(comparison=dict(configurations=[str(base/'none.toml'),str(base/'dual.toml')],labels=['free','dual'],working_total_lengths_km=[100.],workers=1))
            two=compare(spec,workers=1);self.assertEqual(len(two['variants']),2);self.assertEqual(len(two['working_lengths']),2)
            free,dual=two['working_lengths'];self.assertAlmostEqual(dual['rate_ratio_to_first'],dual['key_bps']/free['key_bps'])
            spec['comparison']['configurations'].append(str(base/'classical.toml'));spec['comparison']['labels'].append('classical')
            three=compare(spec,workers=1);self.assertEqual(len(three['variants']),3)
            classical=three['working_lengths'][-1];self.assertTrue(classical['upper_estimate'])
            self.assertIn('CLASSICAL_DETECTION_NOISE',comparison_markdown(three))
            spec['comparison']['labels'].append('extra')
            with self.assertRaises(ValueError):compare(spec,workers=1)

if __name__=='__main__':unittest.main()
