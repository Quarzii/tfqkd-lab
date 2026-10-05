"""Measured inputs above the unchanged spectral and protocol kernels.

Direct log-log interpolation/extrapolation is an explicitly requested numerical
input representation, not a fitted noise mechanism. Detector projections do not
prove bounds for a receiver with two different detectors.
"""
from dataclasses import replace
import copy
import numpy as np
from .keyrates import KeyRateParameters, sns_aopp_rate_per_pulse, cal_rate_per_pulse
from .protocol import duty_cycle, phase_error
from .measured_inputs import normalize_psd, MeasurementRejected


def measured_phase(lab):
    return lab.settings.get('phase', {}).get('mode') == 'measured'


def receiver_inputs(detector):
    channels = detector.get('channels')
    if channels is None:
        return None
    if any(k in detector for k in ('efficiency', 'dark_count_rate_hz')):
        raise ValueError('Use scalar detector characteristics OR two detector channels')
    if not isinstance(channels, list) or len(channels) != 2:
        raise ValueError('detector.channels must contain exactly two detector specifications')
    for entry in channels:
        required={'name','efficiency','dark_count_rate_hz'}
        if not required<=set(entry) or set(entry)-required-{'background_count_rate_hz'}:
            raise ValueError('Each detector requires name, efficiency and dark_count_rate_hz')
        if not isinstance(entry['name'], str) or not entry['name'].strip():
            raise ValueError('Detector name must be nonempty')
        if not np.isfinite(entry['efficiency']) or not 0 < entry['efficiency'] <= 1:
            raise ValueError('Each efficiency must be finite and in (0,1]')
        if not np.isfinite(entry['dark_count_rate_hz']) or entry['dark_count_rate_hz'] < 0:
            raise ValueError('Each dark count rate must be finite and nonnegative [Hz]')
        if 'background_count_rate_hz' in entry and (not np.isfinite(entry['background_count_rate_hz']) or entry['background_count_rate_hz']<0):
            raise ValueError('Explicit background count rate must be finite and nonnegative [Hz]')
    if channels[0]['name'] == channels[1]['name']:
        raise ValueError('Detector names must be distinct')
    return channels


def channel_noise(channel):
    # Zhou2023 Supplementary Note3/TableS3: total noise counts include intrinsic dark + scattered background.
    # QKD.ipynb scalar pdc represents the lumped spurious count contribution; absent background is not a measured zero.
    return channel['dark_count_rate_hz']+channel.get('background_count_rate_hz',0.)


def explicit_losses(line):
    supplied = {'loss_a_db', 'loss_b_db'} & set(line)
    if not supplied:
        return None
    if supplied != {'loss_a_db', 'loss_b_db'}:
        raise ValueError('Both loss_a_db and loss_b_db are required')
    if 'attenuation_db_per_km' in line:
        raise ValueError('Use complete arm losses OR per-km attenuation')
    losses = [line['loss_a_db'], line['loss_b_db']]
    if not np.all(np.isfinite(losses)) or min(losses) < 0:
        raise ValueError('Complete arm losses must be finite and nonnegative [dB]')
    return losses


def effective_loss(lab):
    losses = explicit_losses(lab.settings['line'])
    if losses is None:
        a, b = lab.arms()
        # QKD.ipynb calc_sigma_tau_loss; bertaina2024 Sec.IV: equalize to the lossier arm.
        losses = [a*lab.config['keyrate']['attenuation_db_per_km'],
                  b*lab.config['keyrate']['attenuation_db_per_km']]
    # Same author equalization convention, applied to explicit TOTAL losses with insertions.
    return 2*max(losses)


def protocol_rates(lab, sigma, duty):
    p = KeyRateParameters(**lab.config['keyrate'])
    model = sns_aopp_rate_per_pulse if lab.settings['protocol']['name'] == 'SNS-AOPP' else cal_rate_per_pulse
    loss = effective_loss(lab)
    channels = receiver_inputs(lab.settings['detector'])
    if channels is None:
        # Eq.A1/B1, bertaina2024; QKD.ipynb plot_panel_scenarios, unchanged protocol function.
        raw = p.clockrate_hz*duty*model(loss, sigma, p)
        return dict(raw_key_bps=raw, key_bps=np.maximum(0, raw))
    raw = []
    for channel in channels:
        scalar = replace(p, detector_efficiency=channel['efficiency'],
                         detector_dark_count_rate_hz=channel_noise(channel))
        # QKD.ipynb scalar receiver: hypothetical two identical detectors with this channel's characteristics.
        raw.append(p.clockrate_hz*duty*model(loss, sigma, scalar))
    raw = np.stack(raw)
    clipped = np.maximum(0, raw)  # QKD.ipynb plotting convention for a nonpositive secure-rate bound.
    # Engineering display choice: the smaller projection is a reference result, NOT a proven physical lower bound.
    return dict(raw_key_bps=np.min(raw, axis=0), key_bps=np.min(clipped, axis=0),
                projection_0_raw_bps=raw[0], projection_1_raw_bps=raw[1],
                projection_0_bps=clipped[0], projection_1_bps=clipped[1],
                projection_min_bps=np.min(clipped, axis=0), projection_max_bps=np.max(clipped, axis=0))


def receiver_report(lab, result):
    channels = receiver_inputs(lab.settings['detector'])
    losses = explicit_losses(lab.settings['line'])
    report = dict(complete_arm_losses_db=losses, effective_equalized_loss_db=effective_loss(lab),
                  loss_convention='Bertaina Sec.IV / author calc_sigma_tau_loss: attenuate the stronger arm; effective loss = 2*max(arm losses). This models equalization, not arbitrary asymmetric encoding.')
    if channels:
        report.update(channels=[dict(c, protocol_noise_count_rate_hz=channel_noise(c),key_bps=result[f'projection_{i}_bps'],
                                     raw_key_bps=result[f'projection_{i}_raw_bps']) for i, c in enumerate(channels)],
                      projection_range_bps=[result['projection_min_bps'], result['projection_max_bps']],
                      reduction='No detector averaging. Two scalar projections, each assuming both ports have the specified detector characteristics. The smaller projection is displayed as the reference; the interval is not a certified bound for two unequal detectors.')
    else:
        report.update(channels=None, projection_range_bps=None, reduction='Explicit scalar receiver input.')
    return report


def direct_psd(spec, node, base_dir, grid):
    f, psd, metadata = normalize_psd(spec, base_dir)
    if len(f) < 2:
        raise MeasurementRejected('Direct log-log interpolation needs at least two distinct positive PSD points')
    policy = spec.get('extrapolation', 'reject')
    if policy not in ('reject', 'power-law'):
        raise MeasurementRejected('extrapolation must be reject or explicitly chosen power-law')
    if node == 'line':
        length = spec.get('measurement_length_km')
        if length is None or not np.isfinite(length) or length <= 0:
            raise MeasurementRejected('Direct arm PSD needs measurement_length_km of the single-pass measured path')
    outside = bool(grid['f_min_hz'] < f[0] or grid['f_max_hz'] > f[-1])
    if outside and policy == 'reject':
        raise MeasurementRejected(f'Calculation grid {grid["f_min_hz"]}..{grid["f_max_hz"]} Hz exceeds measured band {f[0]}..{f[-1]} Hz. Narrow the grid or explicitly enable power-law extrapolation.')
    # User-requested log-log interpolation; endpoint slopes are numerical secants, not fitted physical coefficients.
    slopes = [float((np.log(psd[1])-np.log(psd[0]))/(np.log(f[1])-np.log(f[0]))),
              float((np.log(psd[-1])-np.log(psd[-2]))/(np.log(f[-1])-np.log(f[-2])))]
    metadata.update(mode='direct', measured_band_hz=[float(f[0]), float(f[-1])],
                    extrapolation=policy, extrapolation_used=outside, endpoint_log_slopes=slopes,
                    interpolation='piecewise linear log(PSD) versus log(f); no F1/Eq.6 fit',
                    interpretation=('measured laser output PSD; no F2 applied' if node == 'laser' else
                                    'single-arm fiber contribution for the selected compensation state, excluding laser and central detection noise; no other compensation state inferred'))
    return dict(frequency_hz=f, phase_psd=psd, extrapolation=policy,
                measurement_length_km=spec.get('measurement_length_km')), metadata


def interpolate_psd(data, frequency):
    frequency = np.asarray(frequency, dtype=float)
    f, psd = data['frequency_hz'], data['phase_psd']
    outside = (frequency < f[0]) | (frequency > f[-1])
    if np.any(outside) and data['extrapolation'] == 'reject':
        raise MeasurementRejected('A refined calculation grid exceeds the measured PSD band; extrapolation was not authorized')
    # User-requested numerical interpolation/extrapolation of measured samples in log-log coordinates.
    x, y, q = np.log(f), np.log(psd), np.log(frequency)
    result = np.interp(q, x, y)
    for mask, start, end in [(q < x[0], 0, 1), (q > x[-1], -2, -1)]:
        if np.any(mask):
            slope = (y[end]-y[start])/(x[end]-x[start])
            result = np.where(mask, y[start]+slope*(q-x[start]), result)
    output = np.exp(result)
    if np.any(~np.isfinite(output)) or np.any(output <= 0):
        raise MeasurementRejected('Direct PSD interpolation/extrapolation overflowed or underflowed; change the explicit frequency band')
    return output


def mixed_components(f, lab, sc):
    from .spectra import free_laser, stabilized_laser, free_fiber, stabilized_fiber, detection
    from .transfers import common_laser_power
    a, b = lab.arms();p = lab.config['physics'];direct = lab.direct_spectra
    laser = (interpolate_psd(direct['laser'], f) if 'laser' in direct else
             (stabilized_laser if sc['cavity'] else free_laser)(f, p))
    if 'line' in direct:
        measured = interpolate_psd(direct['line'], f)
        # Eq.6/8 length dependence, bertaina2024: explicit spatial scaling assumption for uploaded single-arm PSD.
        fiber = measured*(a+b)/direct['line']['measurement_length_km']
    else:
        fn = stabilized_fiber if sc['stabilized'] else free_fiber
        fiber = fn(f, a, p)+fn(f, b, p)
    # Eq.5/7, bertaina2024: insert uploaded arm/laser terms without imposing their parametric shapes.
    laser = laser*common_laser_power(f, a-b, p) if sc['common'] else 2*laser
    if sc['common']:fiber = p['K']*fiber
    det = detection(f, p) if sc['stabilized'] else np.zeros_like(f)
    # Appendix G, bertaina2024: central detection contribution added once.
    return dict(laser=laser, fiber=fiber, detection=det, total=laser+fiber+det)


def phase_result(lab):
    phase = lab.settings['phase'];sigma = phase['sigma_phi_rad'];tau = phase['tau_s']
    # Eq.1, bertaina2024: small-phase error is reported explicitly; protocol retains the unchanged author Gaussian model.
    variance = sigma**2
    d = float(duty_cycle(tau, phase['tau_ps_s']))
    rates = {k:float(v) for k,v in protocol_rates(lab, sigma, d).items()}
    line = lab.settings['line'];geometry = ('arm_a_km' in line or 'length_km' in line)
    arms = list(lab.arms()) if geometry else None
    result = dict(**rates, tau_q_s=tau, variance_rad2=variance, sigma_phi_rad=sigma,
                  e_phi=float(phase_error(variance)), duty=d, status='user-specified measured window; no threshold solve',
                  loss_db=effective_loss(lab), physical_length_km=sum(arms) if arms else None,
                  arm_lengths_km=arms, protocol=lab.settings['protocol']['name'], compensation='measured residual phase',
                  phase_mode='measured', loop=None, classical_diagnostic=None,
                  input_case_calculations=1, spectral_calculations=0,
                  provenance=lab.provenance, warnings=lab.warnings, missing_amplitudes=[],
                  rate_estimate=('upper model estimate: unknown detector amplitude / intrinsic error / fEC replaced by ideal values' if lab.settings.get('_ideal_protocol_fields') or lab.missing_amplitudes else
                                 'estimate conditional on measured residual phase; not a noise forecast'),
                  phase_measurement=copy.deepcopy(phase))
    result['missing_amplitudes']=list(lab.missing_amplitudes)
    result['receiver'] = receiver_report(lab, result)
    return result
