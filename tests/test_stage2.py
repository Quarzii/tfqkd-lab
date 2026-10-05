"""Meaningful stage-2 boundary and independent quadrature checks."""
import unittest
import numpy as np
from tfqkd.config import load
from tfqkd.classical import remote_power_ratio,high_gain_ratio,low_frequency_ratio,delay_and_boundary,integrator_controller,remote_integrator_ratio,dimensionless_poles,phase_resonance_frequency
from tfqkd.linewidth import numerical_linewidth,white_linewidth

class StageTwoTests(unittest.TestCase):
    def test_white_linewidth_from_fourier_integral(self):
        c=load()['linewidth']
        for h0 in c['white_levels_hz2_per_hz']:
            self.assertAlmostEqual(numerical_linewidth(h0,c['quadrature_tolerance'])/white_linewidth(h0),1)

    def test_missing_loop_and_invalid_gain(self):
        c=load(); length=c['classical']['length_km'];tau,_,b=delay_and_boundary(length,c['physics'])
        with self.assertRaises(ValueError):
            integrator_controller(0)
        with self.assertRaises(ValueError):
            remote_power_ratio(np.array([b]),tau,None,c['classical']['spatial_nodes'])

    def test_no_gain_remote_is_free_noise(self):
        c=load();tau,_,b=delay_and_boundary(c['classical']['length_km'],c['physics'])
        # Eqs. A4,A8, williams2008: disabled controller leaves original one-way spectrum.
        ratio=remote_power_ratio(np.array([b]),tau,lambda s: np.zeros_like(s),c['classical']['spatial_nodes'])
        np.testing.assert_allclose(ratio,1,rtol=c['validation']['reference_rtol'])

    def test_zero_delay_ideal_gain(self):
        c=load()
        # Eq. A9, williams2008: high-gain remote residual vanishes for zero delay.
        self.assertEqual(float(high_gain_ratio(np.array([c['physics']['fc1_hz']]),0,c['classical']['spatial_nodes'])[0]),0)


    def test_analytic_spatial_integral(self):
        c=load(); tau,_,b=delay_and_boundary(c['classical']['length_km'],c['physics'])
        f=np.geomspace(b/100,b*10,100)
        for g in c['classical']['comparison_g_per_s']:
            expected=remote_power_ratio(f,tau,integrator_controller(g),c['classical']['spatial_nodes'])
            actual=remote_integrator_ratio(f,tau,g)
            np.testing.assert_allclose(actual,expected,rtol=c['validation']['reference_rtol'])

    def test_poles_satisfy_characteristic_equation(self):
        c=load(); _,rt,_=delay_and_boundary(c['classical']['length_km'],c['physics'])
        for g in c['classical']['comparison_g_per_s']:
            a=g*rt
            z=dimensionless_poles(a,c['classical']['pole_branches'])
            # Eq. A6 with C=g/s: characteristic equation z+a+a*exp(-z)=0.
            np.testing.assert_allclose(np.abs(z+a+a*np.exp(-z))/(1+np.abs(z)),0,atol=c['validation']['reference_rtol'])
            self.assertTrue(np.all(z.real<0))

    def test_resonance_scale_inverse_length(self):
        c=load();products=[]
        for length in c['classical']['scaling_lengths_km']:
            _,rt,_=delay_and_boundary(length,c['physics'])
            products.append(length*phase_resonance_frequency(rt))
        np.testing.assert_allclose(products,products[0],rtol=c['validation']['reference_rtol'])

    def test_large_gain_remote_limit_is_nonzero(self):
        c=load();p=c['classical']
        tau,rt,_=delay_and_boundary(p['length_km'],c['physics'])
        # Eqs. A8-A9, williams2008: compare at a fixed low frequency away from resonances.
        f=np.array([phase_resonance_frequency(rt)*p['low_frequency_fraction']])
        expected=high_gain_ratio(f,tau,p['spatial_nodes'])
        actual=remote_integrator_ratio(f,tau,p['high_gain_dimensionless'][-1]/rt)
        self.assertGreater(float(actual[0]),0)
        # Numerical tolerance for this configured asymptotic probe; no fitted constants.
        np.testing.assert_allclose(actual,expected,rtol=p['quadrature_rtol'],atol=0)

if __name__=='__main__':
    unittest.main()
