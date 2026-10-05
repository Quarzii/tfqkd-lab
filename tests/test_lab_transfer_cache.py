"""Cache scheduling must preserve transfer values and all parameter identities."""
import unittest
from unittest.mock import patch
import numpy as np
from tfqkd.lab_transfer_cache import cached_ratio, clear, search_ratio
from tfqkd.classical import remote_controller_ratio

class TransferCacheTests(unittest.TestCase):
    def setUp(self):clear()

    def test_hoisted_search_is_bitwise_equal_to_original_core(self):
        from tfqkd.classical import actuator_stability
        f=np.geomspace(1e-4,1e7,4097)
        for tau in (0.,1e-6,.001,.01):
            for omega in (2*np.pi*100,2*np.pi*100000):
                critical=actuator_stability(2*tau,omega)['g_crit_per_s'] if tau else omega
                gains=critical*np.array([0.,.001,.25,.5,.999999])
                controller=lambda s:gains[:,None]/(s*(1+s/omega)) # Eq.A6 + adopted actuator pole.
                np.testing.assert_array_equal(search_ratio(f,tau,controller),
                                              remote_controller_ratio(f,tau,controller))

    def test_exact_reuse_and_each_parameter_invalidates(self):
        f=np.geomspace(1.,1e5,129);g=np.array([10.]);tau=.001;omega=1000.
        # Eq.A6, williams2008 + explicitly adopted one-pole engineering controller.
        def response(gain,pole):return lambda s:gain[:,None]/(s*(1+s/pole))
        expected=remote_controller_ratio(f,tau,response(g,omega))
        with patch('tfqkd.lab_transfer_cache.remote_controller_ratio',wraps=remote_controller_ratio) as compute:
            first=cached_ratio(f,tau,response(g,omega),g,omega,1.)
            second=cached_ratio(f.copy(),tau,response(g,omega),g.copy(),omega,1.)
            self.assertIs(first,second);self.assertEqual(compute.call_count,1)
            np.testing.assert_array_equal(first,expected)
            with self.assertRaises(ValueError):first[0,0]=0.
            for frequency,delay,gain,pole in [(f*2,tau,g,omega),(f,tau*2,g,omega),(f,tau,g*2,omega),(f,tau,g,omega*2)]:
                cached_ratio(frequency,delay,response(gain,pole),gain,pole,1.)
            self.assertEqual(compute.call_count,5)

    def test_display_keeps_narrow_extrema(self):
        from tfqkd.lab_report import display_indices
        f=np.geomspace(1.,1e6,10000);y=np.ones(len(f));y[1234]=100.;y[4567]=.01
        keep=display_indices(f,y,bins=100)
        self.assertIn(1234,keep);self.assertIn(4567,keep)
        self.assertIn(0,keep);self.assertIn(len(f)-1,keep)
        self.assertLessEqual(len(keep),402)

    def test_zero_ceiling_skips_scan_but_preserves_zero_gain_result(self):
        from tfqkd.config import ROOT
        from tfqkd.lab_inputs import read_toml,resolve
        from tfqkd.lab_engine import calculate_point, evaluate_gains
        raw=read_toml(ROOT/'examples/bertaina2024_table3.toml')
        raw['scheme']['compensation']='classical';raw['line']['length_km']=800.
        lab=resolve(raw)
        with patch('tfqkd.lab_engine._search_values',side_effect=AssertionError('Zero ceiling needs no gain scan')):
            result=calculate_point(lab)
        sampled=evaluate_gains(lab,np.linspace(0,.5*result['loop']['g_crit_per_s'],5))
        self.assertTrue(np.all(sampled['key_bps']==0))
        self.assertEqual(result['loop']['g_per_s'],0.)
        for key in ('key_bps','raw_key_bps','tau_q_s','variance_rad2','duty'):
            self.assertEqual(result[key],sampled[key][0])

    def test_reused_ceiling_matches_direct_selected_calculation(self):
        from tfqkd.config import ROOT
        from tfqkd.lab_inputs import read_toml,resolve
        from tfqkd.lab_engine import calculate
        from tfqkd.lab_run import run_point
        raw=read_toml(ROOT/'examples/laboratory/diode_classical.toml')
        for compensation in ('none','dual','classical'):
            raw['scheme']['compensation']=compensation
            lab=resolve(raw);direct=calculate(lab);reused=run_point(lab)
            for key in ('key_bps','tau_q_s','duty','variance_rad2'):self.assertEqual(direct[key],reused[key])

    def test_zero_budget_does_not_retain_arrays(self):
        f=np.array([1.,2.]);g=np.array([10.]);omega=1000.
        controller=lambda s:g[:,None]/(s*(1+s/omega)) # Eq.A6 + user-adopted one-pole controller.
        with patch('tfqkd.lab_transfer_cache.remote_controller_ratio',wraps=remote_controller_ratio) as compute:
            cached_ratio(f,.001,controller,g,omega,0.)
            cached_ratio(f,.001,controller,g,omega,0.)
            self.assertEqual(compute.call_count,2)

    def test_search_cache_reuses_exact_array(self):
        f=np.geomspace(1e-4,1e6,129);g=np.array([0.,10.,100.]);omega=1000.
        controller=lambda s:g[:,None]/(s*(1+s/omega)) # Eq.A6 + adopted actuator pole.
        with patch('tfqkd.lab_transfer_cache.search_ratio',wraps=search_ratio) as compute:
            first=cached_ratio(f,.001,controller,g,omega,1.,search=True)
            second=cached_ratio(f,.001,controller,g,omega,1.,search=True)
            self.assertIs(first,second);self.assertEqual(compute.call_count,1)
            np.testing.assert_array_equal(first,remote_controller_ratio(f,.001,controller))

if __name__=='__main__':unittest.main()
