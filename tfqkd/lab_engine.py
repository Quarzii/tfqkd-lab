"""General laboratory geometry and gain selection using unchanged core equations.

Run independently: python -m tfqkd.lab_engine INPUT.toml
Gain search is a numerical scan with bounded local refinement, not a controller
identification or a proof of a global optimum for every possible input.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from .lab_inputs import resolve, validate
from .engine import batch_grid, stability_batch
from .integration import PhaseIntegral, frequency_grid
from .spectra import components, free_fiber, detection
from .lab_transfer_cache import cached_ratio
from .keyrates import KeyRateParameters, sns_aopp_rate_per_pulse, cal_rate_per_pulse
from .protocol import duty_cycle


def scenario(lab):
    a,b = lab.arms()
    # Eq. 5, bertaina2024, delta_L=L_A-L_B.
    return dict(common=lab.settings['scheme']['lasers']=='common',cavity=lab.settings['laser']['model']=='cavity',
                stabilized=lab.settings['scheme']['compensation']=='dual',delta_L_km=a-b)


def stability(lab):
    a,b = lab.arms();p=lab.config['physics'];omega=lab.settings['actuator']['omega_a_rad_s']
    # Eq. A6, williams2008 + user-authorized actuator pole: tau_RT=2*n*L/c, first imaginary-axis crossing.
    rt = 2*p['n']*np.array([a,b])/p['c_km_s']
    gc,x = stability_batch(rt,omega,lab.config['performance'])
    # Eq. A6, williams2008: both independently stabilized arms must be stable at the shared controller gain.
    return dict(g_crit_per_s=float(np.min(gc)),arm_g_crit_per_s=gc.tolist(),
                critical_frequency_hz=(x/(2*np.pi*rt)).tolist())


def _batch(lab,gains,level,peak_grid=True,include_spectra=False,search=False):
    c=lab.config;p=c['physics'];op=c['operation'];a,b=lab.arms();sc=scenario(lab)
    gains=np.atleast_1d(gains).astype(float)
    classical=lab.settings['scheme']['compensation']=='classical' and np.any(gains>0)
    if classical and peak_grid:
        if np.any(gains<=0):
            raise ValueError('Internal batch must contain either positive gains or only open-loop gains')
        omega=lab.settings['actuator']['omega_a_rad_s']
        fa,_=batch_grid(c,np.full(gains.shape,a),gains,omega,level)
        if a==b:
            f=fa
        else:
            fb,_=batch_grid(c,np.full(gains.shape,b),gains,omega,level)
            # Eq. 4 quadrature of the summed A8 arm spectra: resolve peaks of BOTH arms on the union grid.
            f=np.sort(np.concatenate((fa,fb),axis=-1),axis=-1)
    else:
        f=frequency_grid(c,points=(len(frequency_grid(c))-1)*2**level+1)
    if lab.direct_spectra:
        from .lab_observed import mixed_components
        terms=mixed_components(f,lab,sc)
    else:terms=components(f,sc,c,LB_km=b)
    if classical:
        omega=lab.settings['actuator']['omega_a_rad_s']
        # Eq. A6, williams2008, with explicitly user-authorized engineering C=g/[s*(1+s/omega_a)].
        controller=lambda s:gains[:,None]/(s*(1+s/omega))
        # Eq. A8, williams2008 (uniform uncorrelated spatial noise); Eq. 5 or Eq. 7, bertaina2024 arm combination.
        # Applying this remote residual to the Eq. 5 common-laser geometry is an explicit engineering composition,
        # not an independently validated common-laser apparatus model; the original K is preserved.
        ratio_a=cached_ratio(f,p['n']*a/p['c_km_s'],controller,gains,omega,c['performance']['batch_memory_mb'],search=search)
        # Eq.A8, williams2008: equal arms, delay and controller have IDENTICAL ratios; reuse the same calculation.
        ratio_b=ratio_a if a==b else cached_ratio(f,p['n']*b/p['c_km_s'],controller,gains,omega,c['performance']['batch_memory_mb'],search=search)
        fiber=free_fiber(f,a,p)*ratio_a
        fiber+=free_fiber(f,b,p)*ratio_b
        terms['fiber']=(p['K'] if sc['common'] else 1)*fiber
        terms['total']=terms['laser']+terms['fiber']+terms['detection']
    if lab.settings.get('_perfect_fiber',False):
        # User-requested counterfactual ceiling; Eq.5/7, bertaina2024: remove ONLY the fiber summand.
        terms['fiber']=np.zeros_like(terms['fiber'])
        terms['total']=terms['laser']+terms['detection']
    if terms['total'].ndim==1:
        terms={key:np.broadcast_to(value,gains.shape+value.shape) for key,value in terms.items()}
    integral=PhaseIntegral(f,terms['total'])
    # Eq. 4, bertaina2024: solve threshold with the same accumulated-integral interpolation as part A.
    tau,status=integral.operating_time(op['sigma_limit_rad'],op['tau_max_s'])
    variance=integral.variance(tau)
    duty=duty_cycle(tau,op['tau_ps_s'])
    params=KeyRateParameters(**c['keyrate'])
    # QKD.ipynb calc_sigma_tau_loss; bertaina2024 Sec. IV paragraph after Fig. 3:
    # higher-transmission arm receives extra attenuation, so effective total loss is twice the longer-arm loss.
    from .lab_observed import effective_loss,protocol_rates
    loss=effective_loss(lab)
    # Eq. A1 / Eq. B1, bertaina2024; author plot_panel_scenarios: clock times duty times per-pulse bound.
    protocol=protocol_rates(lab,np.sqrt(variance),duty)
    raw_rate=protocol['raw_key_bps']
    # QKD.ipynb plot_panel_scenarios: nonpositive protocol lower bounds provide no positive secret key.
    rate=np.maximum(0,raw_rate)
    result=dict(tau_q_s=tau,variance_rad2=variance,cap_variance_rad2=integral.variance(op['tau_max_s']),
                sigma_phi_rad=np.sqrt(variance),duty=duty,key_bps=rate,raw_key_bps=raw_rate,
                loss_db=np.full(gains.shape,loss),status=status,grid_points=np.full(gains.shape,f.shape[-1]),
                refinement_level=np.full(gains.shape,level))
    result.update({key:value for key,value in protocol.items() if key not in ('key_bps','raw_key_bps')})
    if include_spectra:
        result['_frequency_hz']=f
        result['_components']=terms
    return result


def _converged(lab,gains):
    if lab.settings['scheme']['compensation']!='classical' or np.all(np.asarray(gains)==0):
        return _batch(lab,gains,0)
    previous=_batch(lab,gains,0)
    for level in range(1,lab.config['classical_actuator']['max_refinement_levels']+1):
        current=_batch(lab,gains,level)
        # Numerical convergence criterion inherited unchanged from part A; no physical coefficients fitted.
        error=np.maximum.reduce([np.abs(previous[key]/current[key]-1) for key in ('tau_q_s','variance_rad2','cap_variance_rad2')])
        if np.all(error<=lab.config['classical']['quadrature_rtol']):
            current['convergence_error']=error
            return current
        previous=current
    raise ArithmeticError(f'Laboratory quadrature did not converge: max relative change={np.max(error)}')


def _search_values(lab,gains):
    """Temporary log-grid objective, convergence checked; winner verified on A's pole-resolving grid.

    This changes scheduling of the new numerical search, not the core's fast/reference
    integration. If a log grid cannot converge (e.g. narrow peaks), use the unchanged
    adaptive quadrature for that search batch.
    """
    def log_batch(selected,level):
        # Computational memory scheduling: arrays over gain/frequency, bounded before allocating A8 matrices.
        nodes=(len(frequency_grid(lab.config))-1)*2**level+1
        budget=lab.config['performance']['batch_memory_mb']*1024**2
        chunk=max(1,int(budget/(nodes*160)))
        pieces=[_batch(lab,selected[start:start+chunk],level,peak_grid=False,search=True)
                for start in range(0,len(selected),chunk)]
        return {key:np.concatenate([value[key] for value in pieces]) for key in pieces[0]}
    gains=np.atleast_1d(gains)
    previous=log_batch(gains,0)
    output={key:np.empty_like(value) for key,value in previous.items()}
    active=np.arange(len(gains))
    for level in range(1,lab.config['classical_actuator']['max_refinement_levels']+1):
        current=log_batch(gains[active],level)
        # Numerical SEARCH accuracy: key/time, at the requested rate tolerance. Final output is verified
        # using A's stricter tau/variance/cap criteria and pole grids. Far-above-threshold cap variance
        # is not an optimization objective; no reported cap variance is taken from this temporary mesh.
        error=np.maximum.reduce([np.abs(previous[key]-current[key])/np.maximum(abs(current[key]),np.finfo(float).tiny)
                                 for key in ('tau_q_s','key_bps')])
        accepted=(error<=lab.settings['numerics']['gain_rate_rtol']) & (previous['status']==current['status'])
        for key in output:
            output[key][active[accepted]]=current[key][accepted]
        active=active[~accepted]
        if not len(active):
            return output
        previous={key:value[~accepted] for key,value in current.items()}
    fallback=evaluate_gains(lab,gains[active])
    for key in output:
        output[key][active]=fallback[key]
    return output


def evaluate_gains(lab,gains):
    """Vectorize each gain chunk. Only independent input cases may use processes."""
    validate(lab);gains=np.atleast_1d(gains).astype(float)
    if np.any(~np.isfinite(gains)) or np.any(gains<0):
        raise ValueError('Gain must be finite and nonnegative')
    if 'line' in lab.direct_spectra and np.any(gains>0):
        raise ValueError('An uploaded arm residual already describes its measured compensation state; another controller gain cannot be inferred from it')
    if lab.settings['scheme']['compensation']!='classical' or 'line' in lab.direct_spectra:
        return _converged(lab,gains)
    gc=stability(lab)['g_crit_per_s']
    if np.any(gains>=gc):
        raise ValueError(f'Unstable or marginal gain: max g={gains.max():.12g} >= g_crit={gc:.12g} s^-1. Reduce g or change the actuator/delay.')
    # Numerical scheduling only. Estimate worst starting peak-grid size; same part A memory budget.
    c=lab.config;perf=c['performance'];ac=c['classical_actuator']
    peak=perf['fast_peak_core_points'] if c['grid']['mode']=='fast' else ac['resonance_points']
    tails=perf['fast_peak_tail_points'] if c['grid']['mode']=='fast' else ac['resonance_tail_points']
    nodes=2*len(frequency_grid(c))+2*ac['resonance_branches']*(2*peak+4*tails+1)
    chunk=max(1,int(perf['batch_memory_mb']*1024**2/(nodes*160)))
    groups=[]
    for mask in (gains==0,gains>0):
        indices=np.flatnonzero(mask)
        for start in range(0,len(indices),chunk):
            idx=indices[start:start+chunk]
            groups.append((idx,_converged(lab,gains[idx])))
    keys=set.intersection(*(set(value) for _,value in groups))
    result={key:np.empty(gains.shape,dtype=groups[0][1][key].dtype) for key in keys}
    for idx,value in groups:
        for key in keys:
            result[key][idx]=value[key]
    return result


def _one(values,index):
    return {key:np.asarray(value)[index].item() for key,value in values.items()}


def calculate_point(lab):
    """One input case, optionally optimize g under the user's safety constraint."""
    if not hasattr(lab,'settings'):
        lab=resolve(lab)
    from .lab_observed import measured_phase,phase_result,receiver_report
    if measured_phase(lab):
        validate(lab)
        return phase_result(lab)
    validate(lab);scheme=lab.settings['scheme']['compensation'];evaluations=0;history=[]
    if scheme!='classical' or 'line' in lab.direct_spectra:
        result=_one(evaluate_gains(lab,[0]),0)
        loop=dict(mode=('already measured arm PSD; no loop inference' if 'line' in lab.direct_spectra else 'not applicable'),g_per_s=None,g_crit_per_s=None,stability_margin_fraction=None)
        evaluations=1
    else:
        st=stability(lab);gc=st['g_crit_per_s'];k=lab.settings['loop']['safety_fraction']
        def value(g):
            nonlocal evaluations
            out=_one(evaluate_gains(lab,[g]),0);evaluations+=1
            return out
        if 'g_per_s' in lab.settings['loop']:
            gain=lab.settings['loop']['g_per_s'];result=value(gain);mode='manual'
        else:
            num=lab.settings['numerics'];maximum=k*gc
            ceiling_lab=lab.clone()
            ceiling_lab.settings['scheme']['compensation']='none'
            ceiling_lab.settings['_perfect_fiber']=True
            ceiling=_one(evaluate_gains(ceiling_lab,[0]),0);evaluations+=1
            # The existing perfect-fiber ceiling bounds all nonnegative A8 residuals.
            # A zero ceiling makes every clipped objective zero; preserve the original g=0 tie choice.
            zero_ceiling=ceiling['key_bps']==0
            searched={}
            def search_values(gains):
                nonlocal evaluations
                gains=np.atleast_1d(gains)
                missing=np.unique([g for g in gains if g not in searched])
                if len(missing):
                    fresh=_search_values(lab,missing);evaluations+=len(missing)
                    for index,g in enumerate(missing):searched[g]=_one(fresh,index)
                return {key:np.asarray([searched[g][key] for g in gains]) for key in searched[gains[0]]}
            best_rate=None;search_converged=False
            for refinement in (() if zero_ceiling else range(num['gain_refinements']+1)):
                # Computational scan over the user-authorized stable domain; not a physical range or fitted prior.
                count=(num['gain_samples']-1)*2**refinement+1
                lower=min(lab.config['classical']['g_min_per_s'],maximum)
                gains=np.unique(np.r_[np.linspace(0,maximum,count),np.geomspace(lower,maximum,count)])
                values=search_values(gains)
                index=int(np.argmax(values['key_bps']));candidates=[(float(gains[index]),_one(values,index))]
                # Refine every strict interior sampled maximum, without assuming the curve is unimodal.
                maxima=np.flatnonzero((values['key_bps'][1:-1]>values['key_bps'][:-2]) &
                                     (values['key_bps'][1:-1]>=values['key_bps'][2:]))+1
                for idx in maxima:
                    memo={}
                    def objective(g):
                        nonlocal evaluations
                        memo[g]=_one(search_values([g]),0)
                        return -memo[g]['key_bps']
                    opt=minimize_scalar(objective,bounds=(gains[idx-1],gains[idx+1]),method='bounded',
                                        options={'xatol':num['gain_xatol_fraction']*gc})
                    if not opt.success:
                        raise ArithmeticError('Gain refinement did not converge')
                    candidates.append((float(opt.x),memo[opt.x]))
                gain,result=max(candidates,key=lambda pair:pair[1]['key_bps'])
                history.append(dict(samples=len(gains),g_per_s=gain,key_bps=result['key_bps']))
                if best_rate is not None:
                    change=abs(result['key_bps']-best_rate)/max(abs(result['key_bps']),np.finfo(float).tiny)
                    if change<=num['gain_rate_rtol']:
                        search_converged=True;break
                best_rate=result['key_bps']
            if zero_ceiling:
                gain=0.;result=value(gain);search_converged=True
                history.append(dict(perfect_fiber_ceiling_bps=ceiling['key_bps'],
                    reason='All nonnegative fiber residuals have zero clipped key; original zero-gain tie preserved.'))
            if not search_converged:
                raise ArithmeticError(f'Gain scan optimum has not converged within rate tolerance: {history}')
            verified=result if zero_ceiling else value(gain)
            # Numerical verification against the unchanged part A pole-resolving integration, not a fitted adjustment.
            disagreement=abs(verified['key_bps']-result['key_bps'])/max(abs(verified['key_bps']),np.finfo(float).tiny)
            if disagreement>num['gain_rate_rtol']:
                raise ArithmeticError(f'Search-grid winner disagrees with the core: relative key error={disagreement}')
            result=verified
            history.append(dict(core_winner_verification_relative_key_error=disagreement))
            mode='automatic numerical maximum; scan plus local refinement, resolution checked'
        # Engineering definition of distance to the first stability crossing; not a physical phase margin in degrees.
        loop=dict(**st,mode=mode,g_per_s=float(gain),stability_margin_fraction=float(1-gain/gc),
                  safety_fraction=k,within_automatic_safety_constraint=bool(gain<=k*gc),scan_history=history)
    a,b=lab.arms()
    if scheme=='classical' and 'line' not in lab.direct_spectra:
        from .lab_diagnostics import above_band_diagnostic
        diagnostic=above_band_diagnostic(lab,loop)
    else:
        diagnostic=None
    messages=list(lab.warnings)
    if scheme=='classical':
        from .lab_limits import CLASSICAL_DETECTION_NOISE
        messages.append(CLASSICAL_DETECTION_NOISE)
    return dict(**result,receiver=receiver_report(lab,result),loop=loop,classical_diagnostic=diagnostic,physical_length_km=a+b,arm_lengths_km=[a,b],protocol=lab.settings['protocol']['name'],
                compensation=scheme,input_case_calculations=1,spectral_calculations=evaluations,
                provenance=lab.provenance,warnings=messages,
                assumptions=['shorter arm attenuated to equalize transmissions (bertaina2024 Sec. IV; author calc_sigma_tau_loss)',
                             'for classical: same controller gain and pole on both arms; uniform uncorrelated spatial noise (Williams A8/A11)',
                             'automatic gain: user-requested safety fraction and numerical scan/refinement, no certified global optimization'])


def calculate(lab, *, point=None):
    """Point result plus conditional equipment requirements for unmeasured amplitudes."""
    if not hasattr(lab,'settings'):lab=resolve(lab)
    if point is None:
        result=calculate_point(lab)
    else:
        import copy
        result=copy.deepcopy(point)
        if (result['compensation'] != lab.settings['scheme']['compensation'] or
            result['protocol'] != lab.settings['protocol']['name'] or
            result['arm_lengths_km'] != list(lab.arms())):
            raise ValueError('Reused calculation must represent the selected installation')
        result['input_case_calculations']=0
        result['spectral_calculations']=0
    from .lab_observed import measured_phase
    if lab.settings.get('_ideal_protocol_fields'):
        from .lab_protocol_requirements import protocol_requirements
        requirements = protocol_requirements(lab)
        result['protocol_requirements'] = requirements
        result['input_case_calculations'] += requirements['input_case_calculations']
        result['spectral_calculations'] += requirements['spectral_calculations']
    if measured_phase(lab):
        if lab.missing_amplitudes:
            from .lab_requirements import equipment_requirements
            requirements=equipment_requirements(lab,result);result['equipment_requirements']=requirements
            result['input_case_calculations']+=requirements['input_case_calculations']
        return result
    result['rate_estimate']='upper estimate: unmeasured noise amplitudes set to zero' if lab.missing_amplitudes else 'model point estimate from explicit inputs'
    if lab.settings.get('_ideal_protocol_fields'):
        result['rate_estimate'] = 'upper model estimate: unreported '+', '.join(lab.settings['_ideal_protocol_fields'])+' evaluated at ideal limiting values'
    if lab.settings['scheme']['compensation']=='classical':
        result['rate_estimate']+='; optimistic classical estimate: phase-detection noise omitted'
    result['missing_amplitudes']=list(lab.missing_amplitudes)
    if lab.missing_amplitudes:
        from .lab_requirements import equipment_requirements
        requirements=equipment_requirements(lab,result)
        result['equipment_requirements']=requirements
        result['input_case_calculations']+=requirements['input_case_calculations']
        result['spectral_calculations']+=requirements['spectral_calculations']
    return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path)
    args=parser.parse_args();print(json.dumps(calculate(resolve(args.config)),indent=2))
