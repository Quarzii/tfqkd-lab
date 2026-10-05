"""Reach from an explicitly specified length scan of the unchanged wrapper.

Multiplying R(L) by two does not move its zero. Target-rate reach does move;
the two directions need not give identical distances. No point observation is
silently extrapolated to another length.
"""
import numpy as np
from scipy.optimize import brentq


def _dominant_coordinates(coordinates,monotone):
    def dominates(a,b):
        return (set(a)==set(b) and any(a[k]>b[k] for k in a if k in monotone)
                and all(a[k]>=b[k] if k in monotone else a[k]==b[k] for k in a))
    return [b for b in coordinates if not any(dominates(a,b) for a in coordinates)]


def _length_point(lab,coordinates,imbalance,length,monotone=()):
    from .lab_engine import calculate_point
    probe=lab.clone()
    line=probe.settings['line'];line.pop('arm_a_km',None);line.pop('arm_b_km',None)
    # Eq.5 geometry, bertaina2024: same fixed-imbalance length convention as Part C.
    line.update(length_km=float(length),imbalance_km=imbalance)
    primary=_dominant_coordinates(coordinates,set(monotone))
    points=[]
    def calculate(coordinate):
        candidate=probe.clone()
        for path,value in coordinate.items():candidate.set_parameter(path,value)
        return calculate_point(candidate)
    for coordinate in primary:points.append(calculate(coordinate))
    # User-prescribed B4 monotonic endpoint propagation. Only a positive lower
    # key estimate permits pruning dominated endpoints; signed negative raw
    # rates need not be monotone, so their original envelope is fully evaluated.
    if min(p['raw_key_bps'] for p in points)<=0:
        for coordinate in coordinates:
            if coordinate not in primary:points.append(calculate(coordinate))
    # Existing B4 endpoint/corner envelope; no input distribution is added.
    point=dict(key_bps=min(p['key_bps'] for p in points),raw_key_bps=min(p['raw_key_bps'] for p in points))
    if not np.isfinite(point['raw_key_bps']):raise ArithmeticError('Nonfinite rate encountered during length bracketing; no reach claimed')
    return point,len(points),sum(p['spectral_calculations'] for p in points)


def _crossing_job(job):
    lab,coordinates,imbalance,cache,field,level,a,b,monotone=job
    cache=dict(cache);added={};cases=spectra=0
    def residual(length):
        nonlocal cases,spectra
        if length not in cache:
            point,n,m=_length_point(lab,coordinates,imbalance,length,monotone)
            cache[length]=point;added[length]=point;cases+=n;spectra+=m
        rate=cache[length][field]
        # Numerical conditioning only; the positive-rate root is unchanged.
        return np.log1p(rate/level)-np.log(2.) if level>0 else rate-level
    root=brentq(residual,a,b,xtol=lab.settings['requirements']['root_rtol'],rtol=lab.settings['requirements']['root_rtol'])
    return dict(root_km=float(root),points=added,input_case_calculations=cases,spectral_calculations=spectra)


def reach_report(lab, range_cases=None, *, workers=None, monotone_parameters=()):
    options = lab.settings.get('reach', {})
    target = lab.settings['requirements'].get('target_key_bps')
    result = dict(zero_rate_reach_km=None, target_key_bps=target,
                  input_case_calculations=0, spectral_calculations=0,
                  message='A factor of 2 in key rate corresponds to an undefined distance in reach for this configuration.')
    enabled = options.get('enabled', False)
    if not isinstance(enabled, bool):
        raise ValueError('reach.enabled must be true or false')
    if not enabled:
        result.update(status='disabled', message='Reach scan is off. Set reach.enabled=true to calculate reach where a length dependence is specified.')
        return result
    from .lab_observed import measured_phase
    if measured_phase(lab):
        result['status'] = 'not identified: one measured phase/window does not supply phase versus length'; return result
    if 'loss_a_db' in lab.settings['line']:
        result['status'] = 'not identified: complete arm losses at one point do not supply losses versus length'; return result
    geometry_ranges = [p for p in lab.ranges if p.startswith('line.') and p != 'line.imbalance_km']
    if geometry_ranges:
        result['status'] = 'not identified: fixed working-length scan conflicts with geometric length ranges; specify point geometry';return result
    if lab.ranges and range_cases is None:
        from .lab_uncertainty import calculate_ranges
        envelope = calculate_ranges(lab, workers=1, calculator=calculate_point)
        range_cases = envelope['cases']
        monotone_parameters=envelope['monotone_parameters']
        result['input_case_calculations'] += envelope['input_case_calculations']
        result['spectral_calculations'] += envelope['spectral_calculations']
    coordinates = [r['coordinates'] for r in range_cases] if range_cases else [{}]
    result['rate_curve_basis'] = ('lower sampled input-range envelope; inherits MONOTONICITY_SCOPE/NONMONOTONE_INTERIOR, not certified continuum bounds'
                                  if range_cases else 'explicit point inputs and their labelled ideal limits, if any')
    required = {'minimum_total_km', 'maximum_total_km', 'geometry'}
    bounds = {'minimum_total_km', 'maximum_total_km'} & set(options)
    if bounds and len(bounds) != 2:
        raise ValueError('Supply both reach.minimum_total_km and maximum_total_km, or neither for numerical bracketing')
    if set(options)-required-{'points', 'enabled', 'workers'} or options['geometry'] != 'fixed_imbalance':
        raise ValueError('Reach scan accepts enabled, bounds, points, workers and geometry=fixed_imbalance')
    imbalance = lab.arms()[0]-lab.arms()[1]
    count = options.get('points', lab.settings['reach']['points'])
    if not isinstance(count, int) or count < 3:
        raise ValueError('reach.points must be an integer >= 3')
    cache = {}

    def evaluate(length):
        if length not in cache:
            point,n,m=_length_point(lab,coordinates,imbalance,length,monotone_parameters)
            result['input_case_calculations']+=n;result['spectral_calculations']+=m
            cache[length] = point
        return cache[length]

    if bounds:
        low, high = options['minimum_total_km'], options['maximum_total_km']
        if not np.all(np.isfinite([low, high])) or not abs(imbalance) < low < high:
            raise ValueError('Reach domain must satisfy |imbalance| < minimum_total_km < maximum_total_km')
        domain_source = 'explicit user numerical bounds'
    else:
        # Numerical bracketing from the actual input length; no apparatus length prior or rate law is inserted.
        low = high = sum(lab.arms());floor = abs(imbalance)
        needed = 2*target if target is not None else 0.
        try:
            for _ in range(lab.settings['requirements']['max_bracket_doublings']):
                if evaluate(low)['raw_key_bps'] > needed:break
                low = floor+(low-floor)/2
                if low <= floor:break
            for _ in range(lab.settings['requirements']['max_bracket_doublings']):
                if evaluate(high)['raw_key_bps'] < 0:break
                high = floor+2*(high-floor)
        except (ValueError, ArithmeticError, FloatingPointError) as error:
            result['status'] = 'length search stopped at model/domain limit; no extrapolated reach: '+str(error)
            return result
        if low >= high:
            result['status'] = 'no finite numerical length bracket found';return result
        domain_source = 'automatic numerical bracket from input length; same fixed-imbalance convention as Part C'

    grid = np.linspace(low, high, count)  # Numerical root bracketing; user-specified length domain.
    raw = [evaluate(x)['raw_key_bps'] for x in grid]
    key = [evaluate(x)['key_bps'] for x in grid]
    specifications=[(0.,True)]+([(target,False),(target/2,False),(2*target,False)] if target is not None else [])
    tasks=[];task_keys=[]
    for level,signed in specifications:
        field='raw_key_bps' if signed else 'key_bps'
        values=[evaluate(x)[field]-level for x in grid]
        for a,b,va,vb in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
            if va*vb<0:
                tasks.append((lab,coordinates,imbalance,cache,field,level,a,b,monotone_parameters))
                task_keys.append((level,signed,a,b))
    from .engine import independent_jobs
    worker_count=options['workers'] if workers is None else workers
    if type(worker_count) is not int or worker_count<1:
        raise ValueError('reach.workers must be a positive integer')
    roots_by_bracket={}
    for identifier,answer in zip(task_keys,independent_jobs(_crossing_job,tasks,worker_count)):
        roots_by_bracket[identifier]=answer['root_km'];cache.update(answer['points'])
        result['input_case_calculations']+=answer['input_case_calculations']
        result['spectral_calculations']+=answer['spectral_calculations']
    result['crossing_workers']=worker_count
    result['endpoint_scheduling']='B4 highest-noise endpoint per fixed corner for positive rates only; full signed envelope at nonpositive rates. Inherits the reported monotonicity scope.'

    def crossings(level, signed=False):
        field = 'raw_key_bps' if signed else 'key_bps'
        roots = []
        values = [evaluate(x)[field]-level for x in grid]
        for a, b, va, vb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if va == 0:roots.append(float(a))
            if va*vb < 0:
                roots.append(roots_by_bracket[(level,signed,a,b)])
        if values[-1] == 0:roots.append(float(high))
        roots = sorted(set(roots))
        if max(values) <= 0:
            return dict(threshold_bps=level, roots_km=roots, reach_km=None,
                        status='no rate above the threshold in the explicit domain', multiple_crossings=len(roots)>1)
        status = ('right-censored: rate still meets threshold at explicit upper length bound' if values[-1] > 0 else
                  'crossing found within explicit length domain' if roots else 'no threshold crossing in explicit domain')
        return dict(threshold_bps=level, roots_km=roots, reach_km=(max(roots) if roots and values[-1] <= 0 else None),
                    status=status, multiple_crossings=len(roots)>1)

    zero = crossings(0., signed=True)
    result.update(status='evaluated on source-model length scan', geometry=options['geometry'],
                  domain_total_km=[low, high], domain_source=domain_source, zero_rate=zero, zero_rate_reach_km=zero['reach_km'],
                  sampled_curve=dict(length_total_km=grid.tolist(), key_bps=key, raw_key_bps=raw),
                  limitation='Scan brackets sampled crossings; hidden internal crossings are not certified. Same apparatus coefficients and imbalance; automatic g is re-optimized at each length. Uploaded arm length scaling uses the documented Eq.6/8 model assumption.')
    if target is None:
        result['target_status'] = 'supply requirements.target_key_bps to convert rate uncertainty into kilometres'; return result
    levels = dict(nominal=crossings(target), doubled_rate=crossings(target/2), halved_rate=crossings(target*2))
    result['target_reach'] = levels
    lengths = {k:v['reach_km'] for k,v in levels.items()}
    if all(v is not None for v in lengths.values()):
        # User-approved operational definition: invert the actual R(L) at R_target, R_target/2 and 2*R_target.
        up = lengths['doubled_rate']-lengths['nominal'];down = lengths['nominal']-lengths['halved_rate']
        result.update(doubling_shift_km=up, halving_shift_km=down,
            message=f'A factor of 2 in key rate corresponds to {up:.9g} km (doubling) / {down:.9g} km (halving) in reach for this configuration.')
    else:
        result['target_status'] = 'one or more target crossings are not identified inside the explicit length domain; no distance extrapolation'
    return result
