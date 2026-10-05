"""B checks: unit conversion, source honesty, geometry, stability and endpoint routing."""
import copy
import unittest
from tests.gates import full_only
import warnings
import numpy as np
from tfqkd.config import load
from tfqkd.lab_inputs import resolve as universal_resolve,reference_input,lorentz_width_to_r2,read_toml
from tfqkd.lab_engine import calculate,evaluate_gains,stability
from tfqkd.lab_uncertainty import monotonicity_audit,calculate_ranges
from tfqkd.integration import frequency_grid,PhaseIntegral
from tfqkd.spectra import components
from tfqkd.engine import converged_batch
from tfqkd.keyrates import KeyRateParameters,rates_for_loss

def resolve(raw=None, **kwargs):
    # Explicit published fixture only in tests; production resolve never uses these values.
    c=load();raw=copy.deepcopy(raw or {})
    laser=raw.get('laser',{})
    if 'preset' in laser:laser['model']='cavity' if laser.pop('preset')=='rio_cavity' else 'free'
    if 'line' in raw and 'length_km' in raw['line'] and not {'arm_a_km','arm_b_km'} & raw['line'].keys():
        raw['line'].setdefault('imbalance_km',.020) # Table I fixture, not a tool default.
    if 'lorentz_width_hz' in laser:
        fixture=read_toml('examples/bertaina2024_table3.toml');fixture['physics'].pop('r2')
        fixture.update({k:v for k,v in raw.items() if k in ('laser','line')})
        for key in raw.get('physics',{}):fixture['physics'][key]=raw['physics'][key]
        return universal_resolve(fixture)
    return reference_input(c,overrides=raw)

class LabTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config=load();cls.audit=monotonicity_audit(cls.config)

    def test_table_i_direct_core(self):
        c=self.config
        for s in c['scenarios']:
            for protocol,key in [('SNS-AOPP','sns_aopp_bps'),('CAL','cal_bps')]:
                a,b=c['operation']['LB_km']+s['delta_L_km'],c['operation']['LB_km']
                lab=resolve(dict(line=dict(arm_a_km=a,arm_b_km=b),laser=dict(preset='rio_cavity' if s['cavity'] else 'rio_free'),
                                 scheme=dict(lasers='common' if s['common'] else 'independent',compensation='dual' if s['stabilized'] else 'none'),
                                 protocol=dict(name=protocol)))
                actual=calculate(lab);loss=2*a*c['keyrate']['attenuation_db_per_km']
                direct=rates_for_loss(c,s,loss)
                for value,reference,tol in [(actual['tau_q_s'],direct['tau_q_s'],c['performance']['regression_time_rtol']),
                                            (actual['variance_rad2'],direct['sigma_phi_rad']**2,c['performance']['regression_variance_rtol']),
                                            (actual['key_bps'],max(0,direct[key]),c['performance']['regression_key_rtol'])]:
                    self.assertLessEqual(abs(value-reference),tol*max(abs(reference),np.finfo(float).tiny))

    def test_white_linewidth_component(self):
        r2=self.config['physics']['r2']
        self.assertAlmostEqual(lorentz_width_to_r2(np.pi*r2),r2)
        lab=resolve(dict(line=dict(length_km=200),laser=dict(preset='rio_free',lorentz_width_hz=np.pi*r2)))
        self.assertEqual(lab.config['physics']['r2'],r2)
        with self.assertRaises(ValueError):resolve(dict(line=dict(length_km=200),laser=dict(linewidth_hz=1000)))
        with self.assertRaises(ValueError):resolve(dict(line=dict(length_km=200),laser=dict(lorentz_width_hz=100),physics=dict(r2=300)))

    def test_geometry_swap_and_loss_equalization(self):
        lab=resolve(dict(line=dict(length_km=200,imbalance_km=2.5)))
        swapped=resolve(dict(line=dict(length_km=200,imbalance_km=-2.5)))
        a=calculate(lab);b=calculate(swapped)
        self.assertAlmostEqual(a['key_bps']/b['key_bps'],1,places=12)
        self.assertEqual(a['physical_length_km'],200)
        self.assertEqual(a['loss_db'],202.5*self.config['keyrate']['attenuation_db_per_km'])
        with self.assertRaises(ValueError):resolve(dict(line=dict(length_km=2,imbalance_km=2.5)))
        with self.assertRaises(ValueError):resolve(dict(line=dict(length_km=200,arm_a_km=100,arm_b_km=100)))

    def test_published_detector_table_ii(self):
        raw=read_toml('examples/bertaina2024_table3.toml')
        raw['keyrate'].pop('detector_efficiency');raw['keyrate'].pop('detector_dark_count_rate_hz')
        raw['detector']=dict(efficiency=.25,dark_count_rate_hz=50) # Table II, bertaina2024 SPAD example.
        lab=universal_resolve(raw)
        self.assertEqual(lab.config['keyrate']['detector_efficiency'],.25)
        self.assertEqual(lab.config['keyrate']['detector_dark_count_rate_hz'],50)
        raw['detector']['efficiency']=1.1
        with self.assertRaises(ValueError):universal_resolve(raw)

    def test_missing_actuator_not_filled(self):
        raw=read_toml('examples/bertaina2024_table3.toml');raw.pop('actuator')
        raw['scheme']['compensation']='classical'
        with self.assertRaisesRegex(ValueError,'omega_a'):universal_resolve(raw)
        raw['actuator']=dict(omega_a_rad_s=self.config['classical_actuator']['omega_a_min_rad_s'])
        manual=universal_resolve(raw)
        self.assertEqual(manual.settings['actuator']['omega_a_rad_s'],raw['actuator']['omega_a_rad_s'])

    def test_classical_core_on_balanced_arms(self):
        c=self.config
        lab=resolve(dict(line=dict(length_km=2*c['classical']['length_km'],imbalance_km=0),laser=dict(preset='rio_cavity'),scheme=dict(lasers='independent',compensation='classical')))
        gain=.5*stability(lab)['g_crit_per_s'];omega=lab.settings['actuator']['omega_a_rad_s']
        direct=converged_batch(c,np.array([c['classical']['length_km']]),np.array([gain]),omega)
        actual=evaluate_gains(lab,[gain])
        for key,core in [('tau_q_s','tau_q_s'),('variance_rad2','variance_rad2'),('key_bps','sns_bps')]:
            np.testing.assert_allclose(actual[key],direct[core],rtol=c['performance']['regression_key_rtol'],atol=0)
        with self.assertRaisesRegex(ValueError,'g_crit'):evaluate_gains(lab,[stability(lab)['g_crit_per_s']])

    @full_only
    def test_gain_search_constraints_and_reference(self):
        lab=resolve(dict(line=dict(length_km=100,imbalance_km=0),laser=dict(preset='rio_cavity'),scheme=dict(lasers='independent',compensation='classical')))
        result=calculate(lab);loop=result['loop']
        self.assertLessEqual(loop['g_per_s'],loop['safety_fraction']*loop['g_crit_per_s'])
        reference=lab.clone();reference.config['grid']['mode']='reference'
        reference.settings['loop']['g_per_s']=loop['g_per_s'];ref=calculate(reference)
        for key,tol in [('tau_q_s',1e-4),('variance_rad2',1e-4),('key_bps',1e-3)]:
            self.assertLessEqual(abs(result[key]/ref[key]-1),tol)
        probe=evaluate_gains(lab,[0,.1*loop['g_crit_per_s'],.5*loop['g_crit_per_s']])
        self.assertGreaterEqual(result['key_bps'],np.max(probe['key_bps'])*(1-1e-3))

    def test_monotonicity_and_endpoint_reduction(self):
        self.assertTrue(self.audit['passed'],self.audit['failed_parameters'])
        self.assertEqual(len(self.audit['rows']),112)
        # Computational degenerate interval fixture: published point, no invented uncertainty range.
        ranges={'physics.'+key:dict(minimum=self.config['physics'][key],maximum=self.config['physics'][key],sources=['bertaina2024 Table III'])
                for key in ('r3','r2','C4','C3','C2','l','s0')}
        result=calculate_ranges(resolve(dict(line=dict(length_km=200),ranges=ranges)),audit=self.audit,workers=1)
        self.assertEqual(result['input_case_calculations'],2)
        self.assertEqual(result['output_ranges']['key_bps'][0],result['output_ranges']['key_bps'][1])

    def test_zero_imbalance_extremum_is_included(self):
        # Table I imbalance magnitude; sign change simply exchanges arms. The midpoint has higher key rate.
        lab=resolve(dict(line=dict(length_km=200),ranges={'line.imbalance_km':dict(minimum=-2.5,maximum=2.5,sources=['bertaina2024 Table I, delta_L=2.5 km; arm exchange'])}))
        out=calculate_ranges(lab,audit=self.audit,workers=1)
        middle=lab.clone();middle.set_parameter('line.imbalance_km',0)
        self.assertEqual(calculate(middle)['key_bps'],out['output_ranges']['key_bps'][1])
        self.assertEqual(out['input_case_calculations'],3)
        self.assertFalse(any('NONMONOTONE_INTERIOR' in w for w in out['warnings']))
        lab.ranges['physics.fc1_hz']=dict(minimum=100.,maximum=100.,sources=['bertaina2024 Table III'])
        other=calculate_ranges(lab,audit=self.audit,workers=1)
        self.assertTrue(any('NONMONOTONE_INTERIOR' in w for w in other['warnings']))

    def test_unsourced_range_and_inactive_geometry_rejected(self):
        with self.assertRaisesRegex(ValueError,'sources'):
            resolve(dict(line=dict(length_km=200),ranges={'physics.l':dict(minimum=44,maximum=44)}))
        with self.assertRaises(ValueError):
            resolve(dict(line=dict(arm_a_km=100,arm_b_km=100),ranges={'line.length_km':dict(minimum=200,maximum=200,sources=['test'])}))
        with self.assertRaises(ValueError):resolve(dict(line=dict(length_km=200),operation=dict(LB_km=100)))
        with self.assertRaises(ValueError):
            resolve(dict(line=dict(length_km=200),ranges={'keyrate.stab_overhead_s':dict(minimum=.001,maximum=.001,sources=['bertaina2024 Table II'])}))

if __name__=='__main__':unittest.main()
