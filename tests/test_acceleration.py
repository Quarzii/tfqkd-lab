"""Fast checks of scheduling, caching and inverse quadrature; physical reference tests remain unchanged."""
import copy
import unittest
from unittest.mock import patch
import numpy as np
from tfqkd.config import load
from tfqkd.integration import PhaseIntegral, frequency_grid
from tfqkd.spectra import free_laser, detection, free_fiber
from tfqkd.cache import clear_spectral_cache, cache_info
from tfqkd.engine import stability_batch, resonance_batch, batch_values
from tfqkd.classical import actuator_stability, actuator_resonance_poles, delay_and_boundary
from tfqkd.classical_keyrate import spectrum_at_arms, rates_from_integral
from tfqkd.keyrates import KeyRateParameters, sns_aopp_rate_per_pulse, cal_rate_per_pulse

class AccelerationTests(unittest.TestCase):
    def test_modes(self):
        c=load()
        self.assertEqual(len(frequency_grid(c)),c['grid']['fast_points'])
        self.assertEqual(len(frequency_grid(c,mode='reference')),c['grid']['reference_points'])
        with self.assertRaises(ValueError):frequency_grid(c,mode='other')

    def test_inverse_trapezoid_without_solver(self):
        c=load();f=frequency_grid(c);psd=detection(f,c['physics'])
        integral=PhaseIntegral(f,psd)
        known=c['operation']['comparison_time_s'] if 'comparison_time_s' in c['operation'] else c['validation']['comparison_times_s'][-2]
        limit=np.sqrt(integral.variance(known))
        with patch('scipy.optimize.brentq',side_effect=AssertionError('Repeated root integration forbidden')):
            tau,status=integral.operating_time(limit,c['operation']['tau_max_s'])
        self.assertEqual(status,'threshold');self.assertAlmostEqual(tau/known,1,places=12)

    def test_batched_integrals_and_padding(self):
        c=load();f=frequency_grid(c);p=free_laser(f,c['physics']);scales=np.array([[1.,2.],[3.,4.]])
        grid=np.broadcast_to(f,scales.shape+(len(f),))
        padded_f=np.concatenate((grid,grid[...,-1:]),axis=-1)
        psd=scales[...,None]*p
        padded_p=np.concatenate((psd,psd[...,-1:]),axis=-1)
        batch=PhaseIntegral(padded_f,padded_p);tau,_=batch.operating_time(c['operation']['sigma_limit_rad'],c['operation']['tau_max_s'])
        for index in np.ndindex(scales.shape):
            scalar=PhaseIntegral(f,psd[index]);expected,_=scalar.operating_time(c['operation']['sigma_limit_rad'],c['operation']['tau_max_s'])
            self.assertAlmostEqual(tau[index]/expected,1,places=12)
            self.assertAlmostEqual(batch.variance(tau)[index]/scalar.variance(expected),1,places=12)

    def test_cache_parameter_and_grid_identity(self):
        c=load();f=frequency_grid(c);clear_spectral_cache()
        first=detection(f,c['physics']);second=detection(f.copy(),c['physics'])
        self.assertIs(first,second);self.assertEqual(cache_info()['hits'],1)
        changed=copy.deepcopy(c['physics']);changed['s0']*=2
        np.testing.assert_allclose(detection(f,changed),2*first,rtol=c['validation']['reference_rtol'])
        shifted=f.copy();shifted[-1]*=2
        self.assertIsNot(detection(shifted,c['physics']),first)
        with self.assertRaises(ValueError):first[0]=0

    def test_vector_stability_and_poles(self):
        c=load();lengths=np.array([[10.,114.,275.]])
        rt=2*c['physics']['n']*lengths/c['physics']['c_km_s'];omega=c['classical_actuator']['omega_a_max_rad_s']
        gc,x=stability_batch(rt,omega,c['performance']);gain=c['classical_actuator']['gain_fractions'][-1]*gc
        z,counts=resonance_batch(rt,gain,omega,c['performance']['fast_peak_branches'],c['performance'])
        for index in np.ndindex(rt.shape):
            old=actuator_stability(rt[index],omega)
            self.assertAlmostEqual(gc[index]/old['g_crit_per_s'],1,places=11)
            reference=actuator_resonance_poles(rt[index],gain[index],omega,c['performance']['fast_peak_branches'])
            nearest=np.nanmin(abs(reference[:,None]-z[index][None,:]),axis=1)
            self.assertLess(np.max(nearest/(1+abs(reference))),c['performance']['pole_residual_rtol'])

    def test_vector_physics_against_scalar_same_mesh(self):
        from tfqkd.engine import batch_grid
        c=load();lengths=np.array([[10.,114.]])
        omega=c['classical_actuator']['omega_a_min_rad_s'];rt=2*c['physics']['n']*lengths/c['physics']['c_km_s']
        gc,_=stability_batch(rt,omega,c['performance']);g=c['classical_actuator']['gain_fractions'][1]*gc
        batch=batch_values(c,lengths,g,omega,'classical',0);f,_=batch_grid(c,lengths,g,omega,0)
        for index in np.ndindex(lengths.shape):
            psd=spectrum_at_arms(f[index],lengths[index],c,'classical',g[index],omega)
            # Duplicate quadrature nodes have zero weight; scalar and batch use the SAME A8 PSD.
            scalar=rates_from_integral(PhaseIntegral(f[index],psd),lengths[index],c)
            for key in ('tau_q_s','variance_rad2','sns_raw_bps','cal_raw_bps'):
                self.assertAlmostEqual(batch[key][index]/scalar[key],1,places=11)

    def test_protocol_broadcast(self):
        c=load();p=KeyRateParameters(**c['keyrate']);loss=np.array([[20.,40.],[60.,80.]])
        sigma=np.array([[0.,.1],[.2,.2]])
        for model in (sns_aopp_rate_per_pulse,cal_rate_per_pulse):
            values=model(loss,sigma,p)
            for index in np.ndindex(loss.shape):
                self.assertAlmostEqual(values[index]/model(loss[index],sigma[index],p),1,places=12)

if __name__=='__main__':unittest.main()
