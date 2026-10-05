"""Part C: prescribed factor-two perturbations and exact core variance decomposition.

No new physical model. Direction choices are engineering analysis settings;
unknown amplitudes use their existing conditional requirements, not sensitivity.
Run: python -m tfqkd.lab_analysis INPUT.toml
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from .lab_engine import calculate_point,_batch
from .lab_inputs import resolve,MONOTONE_CANDIDATES
from .integration import PhaseIntegral
from .engine import independent_jobs

PARAMETER_LABELS={
    'physics.r3':('Laser flicker coefficient r3','rad² Hz²'),
    'physics.r2':('Laser white-frequency coefficient r2','rad² Hz'),
    'physics.C4':('Cavity coefficient C4','rad² Hz³'),
    'physics.C3':('Cavity coefficient C3','rad² Hz²'),
    'physics.C2':('Cavity coefficient C2','rad² Hz'),
    'physics.l':('Fiber noise per km l','rad² Hz/km'),
    'physics.s0':('Dual phase-detection floor s0','rad²/Hz'),
    'keyrate.detector_dark_count_rate_hz':('Detector dark counts','Hz'),
    'keyrate.detector_error':('Detector error','fraction'),
    'keyrate.attenuation_db_per_km':('Line attenuation','dB/km'),
    'keyrate.detector_efficiency':('Detector efficiency','fraction'),
    'keyrate.clockrate_hz':('Pulse frequency','Hz'),
    'actuator.omega_a_rad_s':('Actuator pole omega_a','rad/s'),
    'line.imbalance_km':('Signed arm imbalance','km'),
}

def _perturbation(job):
    lab,path,factor,baseline=job
    section,key=path.split('.')
    value=lab.arms()[0]-lab.arms()[1] if path=='line.imbalance_km' else (lab.config[section][key] if section in ('physics','keyrate','operation') else lab.settings[section][key])
    # User-requested factor-two sensitivity: arithmetic perturbation, not a population estimate or uncertainty interval.
    # User-prescribed bounded-efficiency upgrade: halve the shortfall to one, rather than clamp a doubling.
    changed=1-(1-value)/2 if factor=='halve_shortfall' else value*factor
    operation='halve shortfall to unity' if factor=='halve_shortfall' else f'multiply by {factor:g}'
    row=dict(parameter=path,baseline_value=value,changed_value=changed,operation=operation)
    row['label'],row['units']=PARAMETER_LABELS[path]
    row['display_baseline_value']=value;row['display_changed_value']=changed
    if path=='physics.r2' and 'lorentz_width_hz' in lab.settings['laser']:
        row.update(label='Pure Lorentz linewidth',units='Hz',
                   display_baseline_value=lab.settings['laser']['lorentz_width_hz'],
                   display_changed_value=lab.settings['laser']['lorentz_width_hz']*factor)
    if path=='keyrate.clockrate_hz':
        row['note']='Key rate is approximately proportional to pulse frequency; this gain is not due to noise suppression. Dark-count probability per pulse is recalculated by the unchanged protocol.'
    if path in lab.missing_amplitudes:
        return dict(row,status='not evaluated: unmeasured amplitude; see conditional equipment requirement',key_bps=None)
    if value==0:
        return dict(row,status='no factor-two change: explicit zero',key_bps=baseline['key_bps'],gain_bps=0.,gain_fraction=0.)
    probe=lab.clone()
    try:
        if path=='line.imbalance_km':
            # Eq.5 geometry, bertaina2024: halve imbalance at fixed physical sum, for either input representation.
            probe.settings['line']={**lab.settings['line'], 'length_km':sum(lab.arms()),'imbalance_km':changed}
            probe.settings['line'].pop('arm_a_km',None);probe.settings['line'].pop('arm_b_km',None)
        else:probe.set_parameter(path,changed)
    except ValueError as error:
        return dict(row,status='factor-two change outside admissible model domain',reason=str(error),key_bps=None)
    result=calculate_point(probe)
    # User-defined sensitivity gain, not a protocol equation; signed values preserve any worsening.
    delta=result['key_bps']-baseline['key_bps']
    gain=delta/baseline['key_bps'] if baseline['key_bps']>0 else None
    return dict(row,status='evaluated',key_bps=result['key_bps'],gain_bps=delta,gain_fraction=gain,
                tau_q_s=result['tau_q_s'],duty=result['duty'],g_per_s=result['loop']['g_per_s'],
                input_case_calculations=1,spectral_calculations=result['spectral_calculations'])

def sensitivity(lab,baseline,*,workers=None):
    from .lab_observed import measured_phase
    if measured_phase(lab):
        raise ValueError('Measured phase does not identify spectral equipment coefficients for sensitivity ranking')
    directions=lab.settings['sensitivity']['directions']
    jobs=[];excluded=[]
    # These choices implement the user-prescribed what-if calculation, not a guarantee of a realizable hardware upgrade.
    for path,factor in directions.items():
        if path.startswith('keyrate.') and path.split('.')[1] in lab.settings.get('_ideal_protocol_fields', []):
            excluded.append(dict(parameter=path,status='not evaluated: unreported protocol parameter; see conditional target-rate requirement'));continue
        node=('laser' if path in ('physics.r3','physics.r2','physics.C4','physics.C3','physics.C2') else 'line' if path=='physics.l' else None)
        if node in lab.direct_spectra:
            excluded.append(dict(parameter=path,status=f'not evaluated: {node} direct PSD has no fitted model coefficients'));continue
        if path=='keyrate.attenuation_db_per_km' and 'loss_a_db' in lab.settings['line']:
            excluded.append(dict(parameter=path,status='total arm losses supplied; no per-km coefficient inferred'));continue
        if path in ('keyrate.detector_efficiency','keyrate.detector_dark_count_rate_hz') and lab.settings['detector'].get('channels'):
            excluded.append(dict(parameter=path,status='two explicit detector projections; no scalar receiver coefficient averaged'));continue
        if path=='actuator.omega_a_rad_s' and 'line' in lab.direct_spectra:
            excluded.append(dict(parameter=path,status='already measured arm residual; actuator change cannot be inferred'));continue
        section,key=path.split('.')
        container=lab.config[section] if section in ('physics','keyrate','operation') else lab.settings.get(section,{})
        if key not in container and path!='line.imbalance_km':continue
        if path.startswith('physics.C') and lab.settings['laser']['model']!='cavity':
            excluded.append(dict(parameter=path,status='inactive for free-laser model'));continue
        if path=='physics.s0' and lab.settings['scheme']['compensation']!='dual':
            excluded.append(dict(parameter=path,status='inactive outside dual compensation'));continue
        if path=='actuator.omega_a_rad_s' and lab.settings['scheme']['compensation']!='classical':
            excluded.append(dict(parameter=path,status='inactive outside classical compensation'));continue
        jobs.append((lab.clone(),path,factor,baseline))
    count=lab.settings['numerics']['workers'] if workers is None else workers
    rows=independent_jobs(_perturbation,jobs,count)
    ranked=sorted([r for r in rows if r['status']=='evaluated'],key=lambda r:r['gain_bps'],reverse=True)
    for rank,row in enumerate(ranked,start=1):row['rank']=rank
    all_paths={section+'.'+key for section in ('physics','keyrate','operation') for key in lab.config[section]}
    excluded.extend(dict(parameter=path,status='no approved improvement direction; not ranked') for path in sorted(all_paths-set(directions)-{'physics.c_km_s','physics.K','keyrate.stab_overhead_s','keyrate.duty_bb84'}))
    excluded.extend(dict(parameter=node+'.spectrum',status='direct PSD: no coefficient sensitivity or inverse coefficient requirements without a parametric fit') for node in lab.direct_spectra)
    return dict(ranked=ranked,not_ranked=[r for r in rows if r['status']!='evaluated']+excluded,
                unmeasured_parameters=list(lab.missing_amplitudes),
                interpretation='One parameter at a time; signed gains; fixed manual g if supplied, otherwise re-optimize g with the unchanged search. Not an uncertainty interval or a joint upgrade.',
                input_case_calculations=sum(r.get('input_case_calculations',0) for r in rows),
                spectral_calculations=sum(r.get('spectral_calculations',0) for r in rows))

def variance_contributions(lab,result):
    from .lab_observed import measured_phase
    if measured_phase(lab):
        raise ValueError('Measured phase has no spectral component decomposition')
    gain=result['loop']['g_per_s'] or 0.
    values=_batch(lab,[gain],result['refinement_level'],include_spectra=True)
    f=np.asarray(values['_frequency_hz']);f=f[0] if f.ndim>1 else f
    terms={key:np.asarray(value)[0] for key,value in values['_components'].items()}
    tau=result['tau_q_s'];total=result['variance_rad2'];rows=[];cumulative={}
    for name in ('laser','fiber','detection'):
        # Eq.4, bertaina2024: integrate each additive Eq.5/7/Appendix G term at the SAME solved tau_Q.
        integral=PhaseIntegral(f,terms[name]);variance=float(integral.variance(tau))
        # User-requested relative contribution: arithmetic fraction of the summed variance, not an extra noise model.
        fraction=variance/total if total>0 else None
        rows.append(dict(component=name,variance_rad2=variance,fraction=fraction))
        cumulative[name]=integral.cumulative.tolist()
    summed=sum(r['variance_rad2'] for r in rows)
    relative=abs(summed-total)/max(abs(total),np.finfo(float).tiny)
    if relative>lab.config['validation']['reference_rtol']:
        raise ArithmeticError(f'Variance decomposition disagrees with reported core total: {summed} vs {total}, error={relative}')
    return dict(rows=rows,total_variance_rad2=total,sum_variance_rad2=summed,sum_relative_error=relative,
                dominant_component=max(rows,key=lambda r:r['variance_rad2'])['component'] if total>0 else None,
                frequency_hz=f.tolist(),phase_psd_rad2_per_hz={key:value.tolist() for key,value in terms.items()},
                cumulative_variance_rad2=cumulative,lower_hz=1/tau,
                spectral_calculations=1,input_case_calculations=0,
                detection_status=('Classical phase-detection noise is omitted and unknown; zero here is NOT a measured zero.'
                    if lab.settings['scheme']['compensation']=='classical' else 'Appendix G detection is active only for the dual scheme.'),
                source='Eq.4 with separated Eq.5/7 terms and Appendix G detection, bertaina2024; Williams A8 for classical fiber')

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path)
    args=parser.parse_args();lab=resolve(args.config);point=calculate_point(lab)
    print(json.dumps(dict(sensitivity=sensitivity(lab,point),variance=variance_contributions(lab,point)),indent=2))
