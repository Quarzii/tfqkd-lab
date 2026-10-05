"""Endpoint propagation requested in B4; no invented physical preset ranges.

For parameters not verified monotone, a corner envelope is a sampled envelope,
NOT a rigorous bound over a continuum. Interior extrema are a known limitation.
Run: python -m tfqkd.lab_uncertainty INPUT.toml
"""
from __future__ import annotations
import itertools
import json
from pathlib import Path
import numpy as np
from .config import load
from .lab_inputs import resolve, reference_input, MONOTONE_CANDIDATES
from .lab_engine import calculate
from .engine import independent_jobs


def monotonicity_audit(config=None):
    """Table I checks, both protocols, positive key bounds and fixed-cap variance.

    Factor two is a numerical test perturbation, not a confidence interval or an
    equipment uncertainty. The zero C3 test uses a labelled numerical stress input.
    """
    c=load() if config is None else config
    rows=[];failures=set()
    for scenario in c['scenarios']:
        for protocol in ('SNS-AOPP','CAL'):
            base=reference_input(c,scenario,protocol=protocol)
            for path in MONOTONE_CANDIDATES:
                section,key=path.split('.');probe=base.clone()
                initial=base.config[section][key]
                numerical_stress=False
                if initial==0:
                    # Numerical stress fixture only, dimensions of C3; Eq. F3, bertaina2024 supplies the basis.
                    # C4/fc1 defines a reproducible nonzero probe, NOT a published C3 or an uncertainty bound.
                    initial=c['physics']['C4']/c['physics']['fc1_hz']
                    base_case=base.clone();base_case.set_parameter(path,initial)
                    numerical_stress=True
                else:
                    base_case=base
                low=calculate(base_case);probe.set_parameter(path,2*initial);high=calculate(probe)
                variance_change=high['cap_variance_rad2']-low['cap_variance_rad2']
                key_change=high['key_bps']-low['key_bps']
                # Floating-point comparison only: test does not assign physical tolerances or error bars.
                tolerance=c['validation']['reference_rtol']*max(low['key_bps'],high['key_bps'],np.finfo(float).tiny)
                passed=variance_change>=-c['validation']['reference_rtol']*max(low['cap_variance_rad2'],np.finfo(float).tiny) and key_change<=tolerance
                if not passed:
                    failures.add(path)
                rows.append(dict(scenario=scenario['name'],protocol=protocol,parameter=path,
                                 initial_value=initial,changed_value=2*initial,numerical_zero_coefficient_stress=numerical_stress,
                                 low_key_bps=low['key_bps'],high_key_bps=high['key_bps'],key_change_bps=key_change,
                                 cap_variance_change_rad2=variance_change,passed=bool(passed)))
    return dict(scope='Table I point inputs, factor-two numerical stress checks; no global monotonicity proof',
                rows=rows,failed_parameters=sorted(failures),
                monotone_parameters=sorted(set(MONOTONE_CANDIDATES)-failures),passed=not failures)


def _calculate_case(job):
    lab,coordinates,calculator=job
    for path,value in coordinates.items():
        lab.set_parameter(path,value)
    result=calculator(lab)
    return dict(coordinates=coordinates,result=result)


def calculate_ranges(lab,*,audit=None,workers=None,calculator=calculate):
    if not hasattr(lab,'settings'):
        lab=resolve(lab)
    if not lab.ranges:
        result=calculator(lab)
        return dict(point=result,output_ranges={key:[result[key],result[key]] for key in ('key_bps','tau_q_s','duty','variance_rad2')},
                    input_case_calculations=1,spectral_calculations=result['spectral_calculations'],
                    range_method='point input; no uncertainty interval assigned',corner_parameters=[],range_sources={})
    # The published scenarios are an explicit validation fixture, never apparatus defaults for this lab.
    audit=monotonicity_audit() if audit is None else audit
    monotone=set(audit['monotone_parameters'])
    # Eqs. 1,5, didomenico2010 and Eq. F1, bertaina2024: Lorentz width is a positive linear scaling of r2.
    if 'physics.r2' in monotone:
        monotone.add('laser.lorentz_width_hz')
    lowhigh=[key for key in lab.ranges if key in monotone]
    corners=[key for key in lab.ranges if key not in monotone]
    cases=[]
    coordinates_by_path={path:list(dict.fromkeys([lab.ranges[path]['minimum'],lab.ranges[path]['maximum']])) for path in corners}
    zero_imbalance='line.imbalance_km' in corners and lab.ranges['line.imbalance_km']['minimum']<=0<=lab.ranges['line.imbalance_km']['maximum']
    if zero_imbalance:
        # Eq. 5, bertaina2024: the common laser term vanishes at delta_L=0 (known interior extremum).
        coordinates_by_path['line.imbalance_km']=sorted(set(coordinates_by_path['line.imbalance_km'])|{0.})
    for corner in itertools.product(*(coordinates_by_path[path] for path in corners)):
        fixed=dict(zip(corners,corner))
        for side in (('minimum','maximum') if lowhigh else ('minimum',)):
            coordinates=fixed|{path:lab.ranges[path][side] for path in lowhigh}
            cases.append((lab.clone(),coordinates,calculator))
    # Independent endpoint/corner input cases only. The gain scan in each worker is vectorized.
    count=lab.settings['numerics']['workers'] if workers is None else workers
    values=independent_jobs(_calculate_case,cases,count)
    outputs={key:[min(row['result'][key] for row in values),max(row['result'][key] for row in values)]
             for key in ('key_bps','tau_q_s','duty','variance_rad2','cap_variance_rad2')}
    message='monotone endpoint envelope under the verified Table I monotonicity assumption'
    warnings=list(lab.warnings)
    unresolved=[path for path in corners if not (zero_imbalance and path=='line.imbalance_km')]
    if zero_imbalance:
        message='endpoint envelope including the known delta_L=0 extremum; other coefficients use verified Table I monotonicity'
    if unresolved:
        message='sampled corner envelope; NOT guaranteed continuum bounds for nonmonotone inputs'
        warnings.append('NONMONOTONE_INTERIOR: endpoint corners can miss interior extrema; this sampled envelope must not be interpreted as certified uncertainty bounds or a guaranteed lower key rate.')
    return dict(output_ranges=outputs,input_case_calculations=len(values),
                spectral_calculations=sum(row['result']['spectral_calculations'] for row in values),
                cases=values,range_method=message,monotone_parameters=lowhigh,corner_parameters=corners,
                known_interior_points={'line.imbalance_km':[0.]} if zero_imbalance else {},
                unresolved_nonmonotone_parameters=unresolved,
                monotonicity_failures=audit['failed_parameters'],range_sources=lab.ranges,
                sources_status='range source references supplied by user; not automatically independently verified',warnings=warnings)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path);parser.add_argument('--workers',type=int)
    args=parser.parse_args();print(json.dumps(calculate_ranges(resolve(args.config),workers=args.workers),indent=2))
