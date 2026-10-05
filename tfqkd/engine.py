"""Batch execution of the existing equations; numerical choices in [performance].

The g,L arrays broadcast to (N_g,N_L). Working arrays are bounded chunks of
that product with a final frequency axis. Only independent omega_a/scenario
jobs use processes; a vectorized chunk never creates child processes.
"""
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.special import lambertw
from .config import cli_config
from .integration import PhaseIntegral, grid_points
from .classical import remote_controller_ratio
from .spectra import stabilized_laser, free_fiber, stabilized_fiber, detection
from .keyrates import KeyRateParameters, sns_aopp_rate_per_pulse, cal_rate_per_pulse
from .protocol import duty_cycle
from .cache import configure_cache


def stability_batch(rt, omega_a, settings):
    rt = np.asarray(rt, dtype=float)
    if np.any(rt <= 0) or omega_a <= 0:
        raise ValueError('Require positive delay and actuator pole')
    # Eq. A6, williams2008 + authorized actuator pole: x*tan(x/2)=b, b=omega_a*tau_RT.
    b = omega_a * rt
    lo, hi = np.zeros_like(b), np.full_like(b, np.pi)
    for _ in range(settings['stability_iterations']):
        mid = (lo + hi) / 2
        below = mid * np.tan(mid / 2) < b
        lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
    x = (lo + hi) / 2
    # Eq. A6 crossing equation, unchanged from classical.actuator_stability.
    return (x*x + b*b) / (2*b*rt), x


def resonance_batch(rt, gain, omega_a, branches, settings):
    """Same characteristic equation and seeds as the scalar reference; damped vector Newton."""
    rt, gain = np.broadcast_arrays(rt, gain)
    gc, x = stability_batch(rt, omega_a, settings)
    if np.any(~np.isfinite(gain)) or np.any(gain <= 0) or np.any(gain >= gc):
        raise ValueError('Unstable or marginal controller in batch; require 0 < g < g_crit')
    # Eq. A6 + adopted actuator pole: z=s*tau_RT, a=g*tau_RT, b=omega_a*tau_RT.
    a, b = (gain*rt)[..., None], (omega_a*rt)[..., None]
    q = (2*np.arange(branches)+1)*np.pi
    z0 = 1j*q
    # Numerical seeds derived from Eq. A6; the same delayed zeros and low-delay quadratic.
    delayed = z0 - (z0*(1+z0/b)+a*(1+np.exp(-z0))) / (1+a+2*z0/b)
    discriminant = ((1-a)**2 - 8*a/b).astype(complex)
    quadratic = np.concatenate((b*(-(1-a)+np.sqrt(discriminant))/2,
                                 b*(-(1-a)-np.sqrt(discriminant))/2), axis=-1)
    # Ideal-integrator A6 roots used ONLY as initial guesses, as in the original scalar code.
    safe_a = np.minimum(a, 600)
    ideal = lambertw(-safe_a*np.exp(safe_a), np.arange(-branches, branches+1))-safe_a
    ideal = np.where(a <= 600, ideal, np.nan)
    z = np.concatenate((1j*x[..., None], quadratic, delayed, ideal), axis=-1)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        for _ in range(settings['pole_iterations']):
            exp = np.exp(-z)
            value = z*(1+z/b)+a*(1+exp)
            derivative = 1+2*z/b-a*exp
            step = value/derivative
            # Numerical damping only: prevents overflow from rejected trial seeds.
            step /= np.maximum(1, np.abs(step)/settings['pole_step_limit'])
            z -= step
        value = z*(1+z/b)+a*(1+np.exp(-z))
        scale = 1+np.abs(z*(1+z/b))+np.abs(a*(1+np.exp(-z)))
        valid = np.isfinite(z) & (np.abs(value)/scale <= settings['pole_residual_rtol']) & (z.real < 0) & (z.imag > 0)
    # Sort and remove duplicate roots within each batch row. No g/L Python loops.
    z = np.take_along_axis(z, np.argsort(np.where(valid, z.imag, np.inf), axis=-1), axis=-1)
    valid = np.isfinite(z) & (z.real < 0) & (z.imag > 0)
    with np.errstate(invalid='ignore'):
        unique = np.r_[True, np.ones(z.shape[-1]-1, dtype=bool)] & valid
        unique[..., 1:] &= np.abs(np.diff(z, axis=-1)) > 1e-7*(1+np.abs(z[..., 1:]))
    # Root residual must also hold after reordering (invalid positive trial roots are not retained).
    with np.errstate(over='ignore', invalid='ignore'):
        value = z*(1+z/b)+a*(1+np.exp(-z))
        scale = 1+np.abs(z*(1+z/b))+np.abs(a*(1+np.exp(-z)))
        unique &= np.isfinite(value) & (np.abs(value)/scale <= settings['pole_residual_rtol'])
    z = np.take_along_axis(np.where(unique, z, np.nan), np.argsort(~unique, axis=-1, kind='stable'), axis=-1)
    counts = np.sum(unique, axis=-1)
    return z[..., :int(np.max(counts))], counts


def batch_grid(config, lengths, gains, omega_a, level):
    ac, perf = config['classical_actuator'], config['performance']
    mode = config['grid']['mode']
    # Eq. 4, bertaina2024: integrate no lower than 1/tau_max. Grid density is computational.
    lower, upper = 1/config['operation']['tau_max_s'], config['grid']['f_max_hz']
    base = np.geomspace(lower, upper, (grid_points(config)-1)*2**level+1)
    # Appendix A, williams2008: tau_RT is twice the one-way delay.
    rt = 2*config['physics']['n']*lengths/config['physics']['c_km_s']
    branches = (perf['fast_peak_branches'] if mode == 'fast' else ac['resonance_branches'])*2**level
    z, counts = resonance_batch(rt, gains, omega_a, branches, perf)
    # Eq. A6: real and imaginary parts locate peak widths/centers; no alteration of A8 spectrum.
    center, width = z.imag/(2*np.pi*rt[..., None]), -z.real/(2*np.pi*rt[..., None])
    valid = np.isfinite(center) & (center > lower) & (center < upper) & (width > 0) & (width < center)
    core_n = (perf['fast_peak_core_points'] if mode == 'fast' else ac['resonance_points'])*(2**level)
    tail_n = (perf['fast_peak_tail_points'] if mode == 'fast' else ac['resonance_tail_points'])*2**level
    core = ac['resonance_widths']*width
    core_offsets = core[..., None]*np.linspace(-1, 1, core_n)
    tail_end = ac['resonance_tail_fraction']*center
    # Eq. 4 quadrature: resolve exactly the same core and logarithmic wings as the scalar grid.
    with np.errstate(invalid='ignore', divide='ignore'):
        tails = core[..., None]*(tail_end/core)[..., None]**np.linspace(0, 1, tail_n)
    tails = np.where((core < tail_end)[..., None], tails, np.nan)
    offsets = np.concatenate((core_offsets, tails, -tails, np.zeros(core.shape+(1,))), axis=-1)
    local = center[..., None]+offsets
    good = valid[..., None] & (local > lower) & (local < upper) & np.isfinite(local)
    # Padding at the upper endpoint has zero quadrature weight; it does not add spectral power.
    local = np.where(good, local, upper).reshape(lengths.shape+(-1,))
    f = np.sort(np.concatenate((np.broadcast_to(base, lengths.shape+(len(base),)), local), axis=-1), axis=-1)
    return f, counts


def batch_values(config, lengths, gains, omega_a, scheme, level):
    lengths, gains = np.broadcast_arrays(np.asarray(lengths, dtype=float), np.asarray(gains, dtype=float))
    if scheme == 'classical':
        f, counts = batch_grid(config, lengths, gains, omega_a, level)
    elif scheme == 'dual':
        # Eq. 4, bertaina2024: same frequency band and base grid as classical comparison.
        f = np.geomspace(1/config['operation']['tau_max_s'], config['grid']['f_max_hz'], (grid_points(config)-1)*2**level+1)
        counts = np.zeros(lengths.shape, dtype=int)
    else:
        raise ValueError('Unknown compensation scheme')
    physics, operation = config['physics'], config['operation']
    # Eq. 7, bertaina2024: two balanced arms and independent stabilized lasers.
    laser = 2*stabilized_laser(f, physics)
    if scheme == 'dual':
        # Eq. 8 / Appendix G, bertaina2024: suppressed arm terms twice, detection exactly once.
        psd = laser+2*stabilized_fiber(f, lengths, physics)+detection(f, physics)
    else:
        # Eqs. A3,A6,A8, williams2008: authorized actuator controller, evaluated for ALL gains at once.
        tau = (physics['n']*lengths/physics['c_km_s'])[..., None]
        controller = lambda s: gains[..., None]/(s*(1+s/omega_a))
        ratio = remote_controller_ratio(f, tau, controller)
        # Eqs. 6,7, bertaina2024, same free-fiber input used by the pre-A kernel.
        psd = laser+2*free_fiber(f, lengths, physics)*ratio
    integral = PhaseIntegral(f, psd)
    tau, status = integral.operating_time(operation['sigma_limit_rad'], operation['tau_max_s'])
    # Eq. 4, bertaina2024 and unnumbered duty formula / notebook duty_stabilization.
    variance = integral.variance(tau)
    duty = duty_cycle(tau, operation['tau_ps_s'])
    p = KeyRateParameters(**config['keyrate'])
    # Table II / QKD.ipynb calc_sigma_tau_loss: total length is twice the balanced arm length.
    loss = 2*lengths*p.attenuation_db_per_km
    # Appendix A Eq. A1 / Appendix B Eq. B1, bertaina2024: clock and duty multiply the protocol bound.
    sns = p.clockrate_hz*duty*sns_aopp_rate_per_pulse(loss, np.sqrt(variance), p)
    cal = p.clockrate_hz*duty*cal_rate_per_pulse(loss, np.sqrt(variance), p)
    return dict(tau_q_s=tau, operating_status=status, variance_rad2=variance, sigma_phi_rad=np.sqrt(variance),
                duty=duty, loss_db=loss, sns_raw_bps=sns, cal_raw_bps=cal,
                sns_bps=np.maximum(0, sns), cal_bps=np.maximum(0, cal),
                cap_variance_rad2=integral.variance(operation['tau_max_s']),
                resolved_poles=counts, grid_points=np.full(lengths.shape, f.shape[-1]),
                refinement_level=np.full(lengths.shape, level))


def converged_batch(config, lengths, gains, omega_a, scheme='classical'):
    """Refine whole vectorized blocks; no repeated integration in a tau_Q solver."""
    previous = batch_values(config, lengths, gains, omega_a, scheme, 0)
    for level in range(1, config['classical_actuator']['max_refinement_levels']+1):
        current = batch_values(config, lengths, gains, omega_a, scheme, level)
        changes = {name:np.abs(previous[key]/current[key]-1) for name,key in
                   [('time_relative_change','tau_q_s'), ('operating_variance_relative_change','variance_rad2'),
                    ('cap_relative_change','cap_variance_rad2')]}
        error = np.maximum.reduce(list(changes.values()))
        if np.all(error <= config['classical']['quadrature_rtol']):
            current['convergence_error'] = error
            current.update(changes)
            return current
        previous = current
    raise ArithmeticError(f'Batch quadrature did not converge; maximum relative change {np.max(error)}')


def omega_scan(job):
    """One independent omega_a worker. g and L are vectorized, memory-bounded chunks."""
    config, omega_a, lengths, fractions, scheme, gain_mode = job
    perf = config['performance']; configure_cache(perf['cache_entries'])
    lengths = np.asarray(lengths); fractions = np.asarray(fractions)
    # Appendix A A6: critical gains depend on L; the scan samples the requested normalized g.
    rt = 2*config['physics']['n']*lengths/config['physics']['c_km_s']
    gc, _ = stability_batch(rt, omega_a, perf)
    gains = fractions[:, None]*gc[None, :] if gain_mode == 'fraction' else np.broadcast_to(fractions[:, None], (len(fractions), len(lengths)))
    ll = np.broadcast_to(lengths[None, :], gains.shape)
    stable = gains < gc[None, :]
    flat_l, flat_g = ll[stable], gains[stable]
    if not len(flat_l):
        raise ValueError('No stable gain in this job; no spectrum or key rate evaluated')
    # Conservative numerical working-array budget. Chunking changes only scheduling, never the grid.
    ac = config['classical_actuator']
    core = perf['fast_peak_core_points'] if config['grid']['mode']=='fast' else ac['resonance_points']
    tails = perf['fast_peak_tail_points'] if config['grid']['mode']=='fast' else ac['resonance_tail_points']
    branches = perf['fast_peak_branches'] if config['grid']['mode']=='fast' else ac['resonance_branches']
    # Count the roots first, cheaply and for the whole g,L product; grids are still built inside chunks.
    if scheme == 'classical':
        _, pole_counts = resonance_batch(2*config['physics']['n']*flat_l/config['physics']['c_km_s'], flat_g, omega_a, branches, perf)
        roots = 2*int(np.max(pole_counts))
    else:
        roots = 0
    estimated_nodes = 2*grid_points(config)+roots*(2*core+4*tails+1)
    chunk = max(1, int(perf['batch_memory_mb']*1024**2/(estimated_nodes*160)))
    pieces = []
    for start in range(0, len(flat_l), chunk):
        end = start+chunk
        pieces.append(converged_batch(config, flat_l[start:end, None], flat_g[start:end, None], omega_a, scheme))
    values = {}
    for key in pieces[0]:
        data = np.concatenate([piece[key].ravel() for piece in pieces])
        full = np.full(gains.shape, 'unstable' if data.dtype.kind in 'US' else np.nan,
                       dtype='U16' if data.dtype.kind in 'US' else float)
        full[stable] = data
        values[key] = full
    return values | dict(g_per_s=gains, stable=stable)


def scan(config, lengths=None, fractions=None, omegas=None, scheme='classical', workers=None, gain_mode='fraction'):
    ac = config['classical_actuator']
    lengths = np.geomspace(ac['length_min_km'], ac['length_max_km'], ac['length_points']) if lengths is None else lengths
    fractions = ac['gain_fractions'] if fractions is None else fractions
    omegas = np.geomspace(ac['omega_a_min_rad_s'], ac['omega_a_max_rad_s'], ac['omega_a_points']) if omegas is None else omegas
    if gain_mode not in ('fraction', 'fixed'):
        raise ValueError('Gain mode must be fraction or fixed')
    jobs = [(config, float(omega), lengths, fractions, scheme, gain_mode) for omega in omegas]
    workers = config['performance']['workers'] if workers is None else workers
    values = independent_jobs(omega_scan, jobs, workers)
    return {key:np.stack([value[key] for value in values]) for key in values[0]}


def independent_jobs(function, jobs, workers):
    """Only independent scenarios/omega_a/comparisons; never a vectorized frequency axis."""
    if workers == 1 or len(jobs) == 1:
        return list(map(function, jobs))
    else:
        # Only independent omega_a jobs are distributed; each job vectorizes g,L locally.
        with ProcessPoolExecutor(max_workers=workers) as executor:
            return list(executor.map(function, jobs))


if __name__ == '__main__':
    c = cli_config()
    result = scan(c, lengths=[c['classical']['length_km']], omegas=[c['classical_actuator']['omega_a_min_rad_s']])
    print({key:value.tolist() for key,value in result.items()})
