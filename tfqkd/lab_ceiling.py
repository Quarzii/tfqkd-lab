"""User-requested fiber stabilization ceiling, using the unchanged core models.

Detection follows Appendix G: only dual-band compensation adds S_det. Perfect
fiber removes the fiber term from the FREE-fiber baseline, whose S_det is zero.
The user explicitly accepted this convention; no common floor is added.
Run: python -m tfqkd.lab_ceiling INPUT.toml
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from .lab_inputs import resolve, MissingInputs
from .lab_engine import calculate_point

def comparison_ready(lab):
    from .lab_observed import measured_phase
    if measured_phase(lab):
        raise ValueError('A measured phase observation does not identify other compensation schemes or a fiber stabilization ceiling')
    if 'line' in lab.direct_spectra:
        raise ValueError('A direct arm residual describes one measured compensation state; other states and a fiber stabilization ceiling cannot be inferred')
    required=[]
    for key in ('lambda_s_nm','lambda_q_nm','fc2_hz'):
        if key not in lab.config['physics']:required.append('physics.'+key+' (dual comparison)')
    if 'omega_a_rad_s' not in lab.settings['actuator']:required.append('actuator.omega_a_rad_s (classical comparison)')
    if required:raise MissingInputs(required)
    result=lab.clone()
    if 's0' not in result.config['physics']:
        result.config['physics']['s0']=0.
        result.missing_amplitudes+=('physics.s0',)
        result.provenance.append(dict(parameter='physics.s0',value=0.,input='UNMEASURED: optimistic detection contribution in dual comparison'))
        result.warnings.append('DUAL_UPPER_ESTIMATE: detection amplitude s0 is unmeasured; the dual rate uses optimistic zero detection noise.')
    return result

def stabilization_ceiling(lab):
    if not hasattr(lab,'settings'):lab=resolve(lab)
    lab=comparison_ready(lab)
    values={}
    for name,scheme,perfect in [('free','none',False),('perfect','none',True),('classical','classical',False),('dual','dual',False)]:
        probe=lab.clone();probe.settings['scheme']['compensation']=scheme
        probe.settings['_perfect_fiber']=perfect
        # User specification: ceiling comparison uses optimized classical gain, even if the main run requests manual g.
        probe.settings['loop'].pop('g_per_s',None)
        values[name]=calculate_point(probe)
    rates={name:row['key_bps'] for name,row in values.items()}
    # User-requested strict invariant. A violation raises an error; rates are never adjusted to enforce this inequality.
    for name in ('classical','dual'):
        if rates[name]>rates['perfect']:
            raise ArithmeticError(f'Fiber ceiling violated: R_{name}={rates[name]:.12g} > R_perfect={rates["perfect"]:.12g}')
    # User-defined ceiling and realized fractions; undefined/zero-baseline cases are explicit, not arbitrary numbers.
    headroom=rates['perfect']-rates['free']
    factor=rates['perfect']/rates['free'] if rates['free']>0 else ('infinite' if rates['perfect']>0 else None)
    fractions={name:(rates[name]-rates['free'])/headroom if headroom>0 else None for name in ('classical','dual')}
    factor_text=f'{factor:.9g}' if isinstance(factor,float) else str(factor)
    message=(f'Ideal fiber stabilization could increase the key rate by at most a factor of {factor_text} at this length.'
             if factor is not None else 'Ideal fiber stabilization gives no positive key at this length; its gain factor is undefined.')
    warnings=list(lab.warnings)
    from .lab_limits import CLASSICAL_DETECTION_NOISE
    warnings.append(CLASSICAL_DETECTION_NOISE)
    if isinstance(factor,float) and factor<1.1:
        # User-authorized practical decision threshold H<1.1, not a source-derived physical constant.
        warnings.append('Ideal fiber stabilization offers less than a 10% key-rate gain; no fiber compensation scheme can offer a substantial gain under these inputs.')
    return dict(message=message,H=factor,rates_bps={'R_'+name:rate for name,rate in rates.items()},
                realized_ceiling_fraction=fractions,results=values,strict_ceiling_check_passed=True,
                classical_estimate='optimistic: round-trip phase-detection noise omitted',
                classical_realized_fraction_estimate='optimistic; comparison biased in favor of classical compensation',
                detection_convention='S_det only in R_dual (Appendix G); free/classical/perfect baseline has zero S_det',
                unknown_input_interpretation='If amplitudes are unmeasured, these are conditional optimistic comparisons, not an equipment guarantee.',
                input_case_calculations=len(values),spectral_calculations=sum(row['spectral_calculations'] for row in values.values()),warnings=warnings)

def table_i_ceiling_curves(config, arm_lengths_km):
    """Vectorized length scan: no gain optimization needed for H=R_perfect/R_free."""
    from .integration import frequency_grid,PhaseIntegral
    from .spectra import components
    from .keyrates import KeyRateParameters,sns_aopp_rate_per_pulse,cal_rate_per_pulse
    from .protocol import duty_cycle
    lengths=np.asarray(arm_lengths_km,dtype=float)
    if lengths.ndim!=1 or np.any(lengths<=0):raise ValueError('Supply positive shorter-arm lengths')
    f=frequency_grid(config);params=KeyRateParameters(**config['keyrate']);op=config['operation'];rows=[]
    for sc in config['scenarios']:
        # Eq.5/7, bertaina2024: vectorized free-fiber arms, preserving each Table I laser arrangement and imbalance.
        scenario=dict(sc,stabilized=False)
        terms=components(f,scenario,config,LB_km=lengths)
        laser=np.broadcast_to(terms['laser'],terms['total'].shape)
        # User-requested ideal counterfactual: S_F=0, same laser/no detection noise as the free-fiber baseline.
        integrals={'free':PhaseIntegral(f,terms['total']),'perfect':PhaseIntegral(f,laser)}
        # QKD.ipynb calc_sigma_tau_loss: equalize to the longer arm (the defined Table I imbalance is positive).
        loss=2*(lengths+sc['delta_L_km'])*params.attenuation_db_per_km
        for protocol,model in [('SNS-AOPP',sns_aopp_rate_per_pulse),('CAL',cal_rate_per_pulse)]:
            rates={}
            for name,integral in integrals.items():
                tau,_=integral.operating_time(op['sigma_limit_rad'],op['tau_max_s'])
                # Eq.A1/B1, bertaina2024; author plotting code: duty times clock times protocol lower bound.
                rates[name]=np.maximum(0,params.clockrate_hz*duty_cycle(tau,op['tau_ps_s'])*model(loss,np.sqrt(integral.variance(tau)),params))
            if np.any(rates['perfect']<rates['free']):raise ArithmeticError('Vectorized Table I perfect/free ceiling violated')
            H=np.divide(rates['perfect'],rates['free'],out=np.full(lengths.shape,np.nan),where=rates['free']>0)
            rows.append(dict(scenario=sc['name'],protocol=protocol,shorter_arm_length_km=lengths.tolist(),
                             R_free_bps=rates['free'].tolist(),R_perfect_bps=rates['perfect'].tolist(),
                             H=[float(x) if np.isfinite(x) else ('infinite' if p>0 else None) for x,p in zip(H,rates['perfect'])]))
    return rows

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path)
    args=parser.parse_args();print(json.dumps(stabilization_ceiling(resolve(args.config)),indent=2))
