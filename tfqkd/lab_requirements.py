"""User-authorized conditional inverse noise requirements, not equipment defaults.

Run: python -m tfqkd.lab_requirements INPUT.toml
Other unmeasured coefficients are zero for EACH separate requirement. The
requirements cannot be combined into a joint guarantee.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from .lab_inputs import resolve, MONOTONE_CANDIDATES
from .lab_engine import calculate_point

def equipment_requirements(lab,upper=None):
    if not hasattr(lab,'settings'):lab=resolve(lab)
    upper=calculate_point(lab) if upper is None else upper
    baseline=upper['key_bps'];settings=lab.settings['requirements']
    fraction=settings['loss_fraction']
    # User-defined equipment specification, not a protocol equation: R_target=(1-loss_fraction)*R_upper.
    target=(1-fraction)*baseline
    rows=[];count=0;spectra=0
    for path in lab.missing_amplitudes:
        if path not in MONOTONE_CANDIDATES:raise ValueError(f'Inverse search is not authorized for {path}')
        if baseline<=0:
            rows.append(dict(parameter=path,status='undefined: optimistic rate is zero',maximum=None));continue
        # Eq.5, bertaina2024: common-laser phase term is identically zero at delta_L=0.
        # Eq.F3 / Appendix G, bertaina2024: cavity/detection amplitudes are absent from inactive models.
        inactive=(path=='physics.s0' and lab.settings['scheme']['compensation']!='dual')
        inactive|=path in ('physics.C2','physics.C3','physics.C4') and lab.settings['laser']['model']=='free'
        inactive|=path in ('physics.r3','physics.r2','physics.C2','physics.C3','physics.C4') and lab.settings['scheme']['lasers']=='common' and lab.arms()[0]==lab.arms()[1]
        if inactive:
            rows.append(dict(parameter=path,status='unconstrained: coefficient is absent from the selected phase model',maximum=None));continue
        samples={0.:baseline}
        def rate(amplitude):
            nonlocal count,spectra
            if amplitude in samples:return samples[amplitude]
            probe=lab.clone();probe.set_parameter(path,amplitude)
            result=calculate_point(probe);count+=1;spectra+=result['spectral_calculations']
            measured=result['key_bps']
            if not np.isfinite(measured):raise ArithmeticError(f'Nonfinite key model during requirement search for {path}')
            samples[float(amplitude)]=measured
            ordered=sorted(samples)
            values=np.array([samples[x] for x in ordered])
            tolerance=lab.config['performance']['regression_key_rtol']*baseline
            if np.any(np.diff(values)>tolerance):
                raise ArithmeticError(f'Local monotonicity failed for {path}: {list(zip(ordered,values.tolist()))}')
            return measured
        # Numerical bracket seed in the coefficient's own units; it is not a physical estimate or a reported bound.
        high=1.
        physical_end=np.nextafter(lab.config['keyrate']['clockrate_hz'],0.) if path.startswith('keyrate.') else np.inf
        bracketed=False
        for _ in range(settings['max_bracket_doublings']):
            high=min(high,physical_end)
            value=rate(high)
            if value<=target:bracketed=True;break
            if high==physical_end:break
            high*=2
        if not bracketed:
            rows.append(dict(parameter=path,status='no crossing found within numerical search; no finite equipment bound claimed',maximum=None,
                             last_tested_value=max(samples),last_tested_key_bps=value,evaluations=len(samples)-1));continue
        # Resolve sub-unit amplitudes before normalized root solving (e.g. s0 in rad^2/Hz).
        # This is dimensionless numerical bracket contraction, not a physical coefficient prior.
        for _ in range(settings['max_bracket_doublings']):
            lower=high/2
            if rate(lower)>target:break
            high=lower
        else:raise ArithmeticError(f'Could not resolve the inverse bracket for {path}; extend numerical bracket budget')
        # User-authorized monotone one-dimensional inverse search. Normalized bracket avoids dimension-dependent atol.
        normalized=brentq(lambda q:rate(q*high)-target,lower/high,1.,xtol=settings['root_rtol'],rtol=settings['root_rtol'])
        maximum=float(normalized*high);achieved=rate(maximum)
        # User-defined fractional loss; verified by recalculating the unchanged core at the returned coefficient.
        actual_loss=1-achieved/baseline
        if abs(actual_loss-fraction)>lab.config['performance']['regression_key_rtol']:
            raise ArithmeticError(f'Inverse requirement verification failed for {path}: loss={actual_loss}, target={fraction}')
        rows.append(dict(parameter=path,status='conditional single-parameter requirement',maximum=maximum,
                         upper_key_bps=baseline,target_key_bps=target,verified_key_bps=achieved,
                         requested_loss_fraction=fraction,verified_loss_fraction=actual_loss,
                         other_unmeasured_parameters_zero=[p for p in lab.missing_amplitudes if p!=path],evaluations=len(samples)-1,
                         message=f'Your {path} must be below {maximum:.9g} for this link, conditional on the other unmeasured contributions being zero.'))
    return dict(requirements=rows,input_case_calculations=count,spectral_calculations=spectra,
                method='user-authorized monotone one-dimensional search with observed local monotonicity checks',
                interpretation='Separate conditional limits; setting every unknown to its individual maximum does NOT guarantee the requested joint key rate.')

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path)
    args=parser.parse_args();print(json.dumps(equipment_requirements(resolve(args.config)),indent=2))
