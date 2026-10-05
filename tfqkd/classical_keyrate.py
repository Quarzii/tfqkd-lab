"""Adaptive operating time and key rates for Williams A8; user-authorized actuator pole."""
import numpy as np
from .classical import (delay_and_boundary, actuator_controller, actuator_resonance_poles,
                        require_stable_actuator, remote_controller_ratio)
from .config import cli_config
from .integration import PhaseIntegral, grid_points
from .keyrates import KeyRateParameters, sns_aopp_rate_per_pulse, cal_rate_per_pulse
from .protocol import duty_cycle
from .spectra import free_fiber, stabilized_fiber, stabilized_laser, detection


def actuator_grid(config, length_km, g, omega_a_rad_s, level=0):
    settings = config['classical_actuator']
    _, rt, _ = delay_and_boundary(length_km, config['physics'])
    require_stable_actuator(rt, g, omega_a_rad_s)
    # Eq. 4, bertaina2024: the smallest possible lower limit is 1/tau_max, not a fixed working time.
    lower = 1 / config['operation']['tau_max_s']
    upper = config['grid']['f_max_hz']
    count = (grid_points(config)-1)*2**level+1
    frequency = np.geomspace(lower, upper, count)
    # Eq. A6, williams2008 with user pole: characteristic roots locate narrow A8 servo features.
    perf = config['performance']
    fast = config['grid']['mode'] == 'fast'
    branches = (perf['fast_peak_branches'] if fast else settings['resonance_branches']) * 2**level
    poles = actuator_resonance_poles(rt, g, omega_a_rad_s, branches)
    for z in poles:
        # Adopted Eq. A6: pole center and half-width in Hz; these are numerical grid choices.
        center = z.imag / (2 * np.pi * rt)
        width = -z.real / (2 * np.pi * rt)
        if not lower < center < upper or not 0 < width < center:
            continue
        # Eq. 4 quadrature: resolve core and logarithmic wings without changing the PSD.
        core = settings['resonance_widths'] * width
        offsets = np.linspace(-core, core, (perf['fast_peak_core_points'] if fast else settings['resonance_points']) * 2**level)
        tail_end = settings['resonance_tail_fraction'] * center
        if core < tail_end:
            tails = np.geomspace(core, tail_end, (perf['fast_peak_tail_points'] if fast else settings['resonance_tail_points']) * 2**level)
            offsets = np.r_[offsets, tails, -tails]
        local = center + offsets
        frequency = np.r_[frequency, local[(local > lower) & (local < upper)], center]
    return np.unique(frequency), len(poles)


def spectrum_at_arms(frequency, length_km, config, scheme, g=None, omega_a_rad_s=None):
    physics = config['physics']
    # Eq. 7, bertaina2024: balanced arms, two independent cavity-stabilized lasers.
    laser = 2 * stabilized_laser(frequency, physics)
    if scheme == 'dual':
        # Eqs. 7,8 and Appendix G, bertaina2024: length-dependent term twice, detection once.
        return laser + 2 * stabilized_fiber(frequency, length_km, physics) + detection(frequency, physics)
    if scheme != 'classical':
        raise ValueError('Unknown compensation scheme')
    tau, rt, _ = delay_and_boundary(length_km, physics)
    require_stable_actuator(rt, g, omega_a_rad_s)
    # Eqs. A3,A8, williams2008 and Eqs. 6,7, bertaina2024: same free-fiber input and lasers.
    ratio = remote_controller_ratio(frequency, tau, actuator_controller(g, omega_a_rad_s))
    return laser + 2 * free_fiber(frequency, length_km, physics) * ratio


def rates_from_integral(integral, length_km, config):
    params = KeyRateParameters(**config['keyrate'])
    operation = config['operation']
    # Eq. 4 and Sec. V, bertaina2024: find threshold time; 100 ms is only the configured cap.
    tau, status = integral.operating_time(operation['sigma_limit_rad'], operation['tau_max_s'])
    variance = float(integral.variance(tau))
    # Eq. 4, bertaina2024: standard deviation is the square root of integrated variance.
    sigma = float(np.sqrt(variance))
    # Sec. I unnumbered duty formula, bertaina2024 / QKD.ipynb duty_stabilization.
    duty = float(duty_cycle(tau, operation['tau_ps_s']))
    # Table II, bertaina2024 / QKD.ipynb calc_sigma_tau_loss: total distance is two balanced arms.
    loss = float(2 * length_km * params.attenuation_db_per_km)
    # Appendix A Eq. A1 and QKD.ipynb Cell 23: clock and duty multiply SNS-AOPP per-pulse bound.
    sns = float(params.clockrate_hz * duty * sns_aopp_rate_per_pulse(loss, sigma, params))
    # Appendix B Eq. B1 and QKD.ipynb Cell 23: clock and duty multiply CAL per-pulse bound.
    cal = float(params.clockrate_hz * duty * cal_rate_per_pulse(loss, sigma, params))
    if not np.all(np.isfinite([variance, duty, sns, cal])):
        raise ArithmeticError('Nonfinite key-rate bound; outside the numerically validated protocol domain')
    return dict(tau_q_s=tau, operating_status=status, variance_rad2=variance, sigma_phi_rad=sigma,
                duty=duty, loss_db=loss, sns_raw_bps=sns, cal_raw_bps=cal,
                # Appendix A Eq. A1 / Appendix B Eq. B1: a negative lower bound certifies no positive key.
                sns_bps=max(0.0, sns), cal_bps=max(0.0, cal))


def operating_point(config, length_km, scheme, g=None, omega_a_rad_s=None, level=0):
    if config['grid'].get('mode') == 'fast':
        from .engine import batch_values
        values = batch_values(config, np.array([[length_km]]),
                              np.array([[config['classical']['g_min_per_s'] if g is None else g]]),
                              omega_a_rad_s, scheme, level)
        return dict(arm_length_km=float(length_km), scheme=scheme, g_per_s=g, omega_a_rad_s=omega_a_rad_s,
                    **{key:value.item() for key,value in values.items()})
    if scheme == 'classical':
        f, poles = actuator_grid(config, length_km, g, omega_a_rad_s, level)
    else:
        # Eq. 4, bertaina2024: use same band and base numerical density in the dual-band calculation.
        lower = 1 / config['operation']['tau_max_s']
        count = (grid_points(config)-1)*2**level+1
        f = np.geomspace(lower, config['grid']['f_max_hz'], count)
        poles = 0
    psd = spectrum_at_arms(f, length_km, config, scheme, g, omega_a_rad_s)
    integral = PhaseIntegral(f, psd)
    return dict(arm_length_km=float(length_km), scheme=scheme, g_per_s=g, omega_a_rad_s=omega_a_rad_s,
                **rates_from_integral(integral, length_km, config), grid_points=len(f),
                resolved_poles=poles, refinement_level=level,
                # Eq. 4, bertaina2024: diagnostic cap variance, not the actual operating variance.
                cap_variance_rad2=float(integral.variance(config['operation']['tau_max_s'])))


def converged_operating_point(config, length_km, scheme, g=None, omega_a_rad_s=None):
    previous = operating_point(config, length_km, scheme, g, omega_a_rad_s)
    tolerance = config['classical']['quadrature_rtol']
    records = []
    for level in range(1, config['classical_actuator']['max_refinement_levels'] + 1):
        current = operating_point(config, length_km, scheme, g, omega_a_rad_s, level)
        # Eq. 4 quadrature convergence: test working time AND cap integral (narrow peaks may lie below cutoff).
        time_change = abs(previous['tau_q_s'] / current['tau_q_s'] - 1)
        cap_change = abs(previous['cap_variance_rad2'] / current['cap_variance_rad2'] - 1)
        variance_change = abs(previous['variance_rad2'] / current['variance_rad2'] - 1)
        records.append(dict(level=level, time_relative_change=time_change, cap_relative_change=cap_change,
                            operating_variance_relative_change=variance_change))
        if max(time_change, cap_change, variance_change) <= tolerance:
            return current, records
        previous = current
    raise ArithmeticError(f'Grid failed to converge at L={length_km}, g={g}, omega_a={omega_a_rad_s}: {records[-1]}')


if __name__ == '__main__':
    from .classical import actuator_stability
    c = cli_config()
    length = c['classical']['length_km']
    _, rt, _ = delay_and_boundary(length, c['physics'])
    omega_a = c['classical_actuator']['omega_a_min_rad_s']
    gain = c['classical_actuator']['gain_fractions'][0] * actuator_stability(rt, omega_a)['g_crit_per_s']
    print('Engineering assumption: C(s)=g/[s*(1+s/omega_a)]; no apparatus calibration claimed.')
    print(converged_operating_point(c, length, 'classical', gain, omega_a))
    print(converged_operating_point(c, length, 'dual'))
