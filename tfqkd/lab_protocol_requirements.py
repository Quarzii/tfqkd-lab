"""Conditional protocol requirements above the unchanged key-rate kernel.

The target rate is explicit user input. Missing errors are ideal limiting
values for an upper estimate, never inferred properties of the apparatus.
"""
import numpy as np
from scipy.optimize import brentq

# Sec.2.1, Eq.1 discussion, Liu et al., Quantum Frontiers 2,16 (2023), DOI 10.1007/s44214-023-00039-9.
REFERENCE_F_EC = 1.16
REFERENCE_F_EC_SOURCE = 'https://link.springer.com/article/10.1007/s44214-023-00039-9#Sec2'
ERROR_END = np.nextafter(.5, 0.)  # Domain of the existing intrinsic error parameter.


def ideal_protocol_inputs(parameters, settings, provenance, warnings):
    missing = []
    for key, value in [('detector_error', 0.), ('f_error', 1.)]:
        if key not in parameters:
            parameters[key] = value; missing.append(key)
            provenance.append(dict(parameter='keyrate.'+key, value=value,
                input='UNREPORTED: ideal limiting value for an upper estimate, not apparatus default'))
    settings['_ideal_protocol_fields'] = missing
    if missing:
        warnings.append('IDEAL_PROTOCOL_UPPER_ESTIMATE: unreported '+', '.join(missing)+
                        '; ideal limiting values give an upper model estimate, not measured apparatus values.')


def validate_target(settings):
    unknown = set(settings['requirements'])-{'target_key_bps','loss_fraction','root_rtol','max_bracket_doublings'}
    if unknown:
        raise ValueError('Unknown requirements fields: '+', '.join(sorted(unknown)))
    target = settings['requirements'].get('target_key_bps')
    if target is not None and (not np.isfinite(target) or target <= 0):
        raise ValueError('requirements.target_key_bps must be finite and positive [bit/s]')


def protocol_requirements(lab):
    from .lab_engine import calculate_point
    unknown = lab.settings.get('_ideal_protocol_fields', [])
    target = lab.settings['requirements'].get('target_key_bps')
    result = dict(requirements=[], input_case_calculations=0, spectral_calculations=0,
        target_key_bps=target, reference_f_ec=REFERENCE_F_EC, reference_source=REFERENCE_F_EC_SOURCE,
        interpretation='Separate conditional limits, not a joint guarantee; other missing quantities retain their stated optimistic values.')
    if not unknown:
        result['status'] = 'no missing intrinsic error or error-correction factor'; return result
    if target is None:
        result['status'] = 'not computed: supply requirements.target_key_bps [bit/s]; no target is invented'; return result
    result['status'] = 'explicit target evaluated'
    projections = range(2) if lab.settings['detector'].get('channels') else [None]
    tolerance = lab.config['performance']['regression_key_rtol']
    numerical = lab.settings['requirements']
    for missing in unknown:
        for projection in projections:
            probe = lab.clone()
            if missing == 'detector_error':
                probe.set_parameter('keyrate.f_error', REFERENCE_F_EC)
            path = 'keyrate.'+missing
            field = 'key_bps' if projection is None else f'projection_{projection}_bps'
            cache = {}

            def rate(value):
                if value not in cache:
                    probe.set_parameter(path, value)
                    point = calculate_point(probe)
                    result['input_case_calculations'] += 1
                    result['spectral_calculations'] += point['spectral_calculations']
                    cache[value] = point[field]
                    values = [cache[k] for k in sorted(cache)]
                    if np.any(np.diff(values) > tolerance*max(target, max(values))):
                        raise ArithmeticError('Local protocol monotonicity failed for '+path)
                return cache[value]

            low = 0. if missing == 'detector_error' else 1.
            row = dict(parameter=path, detector_projection=projection,
                target_key_bps=target, maximum=None,
                conditional_f_ec=probe.config['keyrate']['f_error'] if missing=='detector_error' else None,
                conditional_e_d=probe.config['keyrate']['detector_error'] if missing=='f_error' else None,
                reference_source=REFERENCE_F_EC_SOURCE if missing=='detector_error' else None)
            if rate(low) < target:
                row.update(status='target unattainable even at the ideal value under these conditions', ideal_key_bps=rate(low))
                result['requirements'].append(row); continue
            high = ERROR_END if missing == 'detector_error' else 2.
            for _ in range(numerical['max_bracket_doublings']):
                if rate(high) <= target:break
                if missing == 'detector_error':break
                high *= 2.  # Numerical bracket expansion, not a physical fEC prior.
            if rate(high) > target:
                row.update(status='no crossing in the tested domain; no finite maximum claimed', last_value=high)
            else:
                root = float(brentq(lambda value:rate(value)-target, low, high,
                                   xtol=numerical['root_rtol'], rtol=numerical['root_rtol']))
                achieved = rate(root)
                relative = abs(achieved/target-1)  # Numerical inverse verification, not a physical equation.
                if relative > tolerance:
                    raise ArithmeticError(f'Protocol requirement verification failed: {path}, {relative}')
                row.update(status='conditional single-parameter requirement', maximum=root,
                           verified_key_bps=achieved, relative_target_error=relative)
            row['evaluations'] = len(cache)
            result['requirements'].append(row)
    return result
