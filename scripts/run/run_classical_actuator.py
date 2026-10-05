"""Williams actuator-pole scan with adaptive tau_Q, duty, and protocol-zero key reach."""

# CLI import bootstrap; no calculation settings are changed.
import sys as _sys
from pathlib import Path as _Path
_PROJECT_ROOT = _Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_PROJECT_ROOT))
ROOT = _PROJECT_ROOT
import csv
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/kvant-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from tfqkd.config import cli_config, ROOT
from tfqkd.classical import (delay_and_boundary, actuator_stability, actuator_pole,
                            actuator_local_sensitivity, remote_controller_ratio, actuator_controller,
                            phase_resonance_frequency, high_gain_ratio)
from tfqkd.classical_keyrate import converged_operating_point
from tfqkd.classical_keyrate import actuator_grid
from tfqkd.keyrates import KeyRateParameters, sns_aopp_rate_per_pulse, cal_rate_per_pulse

OUT = ROOT / 'results'


def csv_write(name, rows):
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (OUT / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    global OUT
    c = cli_config()
    OUT = Path(c.get("_output_dir", OUT))
    p, ac = c['classical'], c['classical_actuator']
    OUT.mkdir(parents=True, exist_ok=True)
    print('Engineering assumption: C(s)=g/[s*(1+s/omega_a)]. omega_a is scanned, not fitted.', flush=True)
    omega_values = np.geomspace(ac['omega_a_min_rad_s'], ac['omega_a_max_rad_s'], ac['omega_a_points'])
    lengths = np.geomspace(ac['length_min_km'], ac['length_max_km'], ac['length_points'])
    reference_length = p['length_km']
    _, reference_rt, _ = delay_and_boundary(reference_length, c['physics'])
    stability_rows, gain_rows, crossing_checks = [], [], []
    for omega_a in omega_values:
        for length in np.unique(np.r_[lengths, reference_length, p['scaling_lengths_km']]):
            _, rt, _ = delay_and_boundary(float(length), c['physics'])
            st = actuator_stability(rt, float(omega_a))
            stability_rows.append(dict(arm_length_km=float(length), omega_a_rad_s=float(omega_a), **st))
        st = actuator_stability(reference_rt, float(omega_a))
        # Eq. A6 with adopted actuator pole: inspect characteristic roots on both sides of marginal gain.
        for fraction in (1 - ac['root_crossing_fraction'], 1 + ac['root_crossing_fraction']):
            a = fraction * st['g_crit_per_s'] * reference_rt
            z = actuator_pole(a, st['b'], 1j * st['x'])
            crossing_checks.append(dict(omega_a_rad_s=float(omega_a), g_fraction=fraction,
                                       root_real_per_s=float(z.real / reference_rt),
                                       root_imag_hz=float(z.imag / (2 * np.pi * reference_rt))))
        gains = np.geomspace(p['g_min_per_s'], ac['gain_ceiling_fraction'] * st['g_crit_per_s'], ac['gain_scan_points'])
        for g in gains:
            gain_rows.append(dict(omega_a_rad_s=float(omega_a), g_per_s=float(g), g_fraction=float(g / st['g_crit_per_s']),
                                  g_crit_per_s=st['g_crit_per_s'], critical_frequency_hz=st['critical_frequency_hz']))
    csv_write('actuator_stability.csv', stability_rows)
    csv_write('actuator_gain_scan.csv', gain_rows)
    convergence, cache = [], {}

    # Stage A scheduling only: independent omega_a/study jobs, vectorized g,L within each job.
    # The equations, curve samples, root criterion and all later diagnostics are unchanged.
    from tfqkd.engine import independent_jobs, omega_scan
    jobs = [(c, float(wa), lengths, ac['gain_fractions'], 'classical', 'fraction') for wa in omega_values]
    jobs += [(c, float(wa), lengths, p['comparison_g_per_s'], 'classical', 'fixed') for wa in omega_values]
    jobs += [(c, float(omega_values[0]), lengths, [ac['gain_fractions'][0]], 'dual', 'fraction')]
    print('Computing vectorized curve jobs:', len(jobs), flush=True)
    batches = independent_jobs(omega_scan, jobs, c['performance']['workers'])
    check_fields = ('time_relative_change', 'cap_relative_change', 'operating_variance_relative_change')
    for job, values in zip(jobs, batches):
        _, wa, sample_lengths, settings, scheme, mode = job
        for i, setting in enumerate(settings):
            for j, length in enumerate(sample_lengths):
                if not values['stable'][i,j]:
                    continue
                # Serialize already computed arrays; these loops perform no spectra or integration.
                _, rt, _ = delay_and_boundary(float(length), c['physics'])
                g = float(setting * actuator_stability(rt, wa)['g_crit_per_s'] if mode == 'fraction' else setting)
                if scheme == 'dual':
                    g, row_wa = None, None
                else:
                    row_wa = wa
                row = {key:value[i,j].item() for key,value in values.items()
                       if key not in (*check_fields, 'convergence_error', 'stable', 'g_per_s')}
                row.update(arm_length_km=float(length), scheme=scheme, g_per_s=g, omega_a_rad_s=row_wa)
                cache[(float(length), scheme, g, row_wa)] = row
                convergence.append(dict(arm_length_km=float(length), scheme=scheme, g_per_s=g,
                                        omega_a_rad_s=row_wa, level=int(row['refinement_level']),
                                        **{key:float(values[key][i,j]) for key in check_fields}))
    print('Vectorized curve points cached:',len(cache),flush=True)

    def evaluate(length, scheme, g=None, omega_a=None):
        key = (float(length), scheme, g, omega_a)
        if key not in cache:
            row, checks = converged_operating_point(c, float(length), scheme, g, omega_a)
            for check in checks:
                convergence.append(dict(arm_length_km=float(length), scheme=scheme, g_per_s=g,
                                        omega_a_rad_s=omega_a, **check))
            cache[key] = row
        return cache[key]

    def root_kwargs():
        return dict(xtol=ac['root_xtol_km'], rtol=ac['root_rtol'])

    def roots_for_curve(rows, callback):
        result = {}
        for protocol in ('sns', 'cal'):
            field = protocol + '_raw_bps'
            roots = []
            # Appendix A Eq. A1 / Appendix B Eq. B1: solve signed key lower bound, not clipped zero plateau.
            for left, right in zip(rows[:-1], rows[1:]):
                if left.get(field) is None or right.get(field) is None:
                    continue
                if left[field] > 0 >= right[field]:
                    root = brentq(lambda length: callback(length)[field], left['arm_length_km'], right['arm_length_km'], **root_kwargs())
                    roots.append(float(root))
            result[protocol + '_zero_arm_km'] = roots[0] if roots else None
            result[protocol + '_crossings_arm_km'] = roots
        return result

    params = KeyRateParameters(**c['keyrate'])
    phase_limit_key_zeros = {}
    for protocol, model in [('sns', sns_aopp_rate_per_pulse), ('cal', cal_rate_per_pulse)]:
        # Appendix A Eq. A1 / Appendix B Eq. B1: duty>0 cannot change a key-bound sign at fixed sigma.
        loss_zero = brentq(lambda loss: model(loss, c['operation']['sigma_limit_rad'], params),
                           ac['key_zero_loss_min_db'], ac['key_zero_loss_max_db'], rtol=ac['root_rtol'])
        # Table II, bertaina2024: convert total channel loss to balanced arm length.
        phase_limit_key_zeros[protocol] = float(loss_zero / (2 * params.attenuation_db_per_km))

    dual_rows = [evaluate(float(length), 'dual') for length in lengths]
    dual_reach = roots_for_curve(dual_rows, lambda length: evaluate(length, 'dual'))
    reference_points = [evaluate(reference_length, 'dual')]
    comparisons, reaches = [], []
    for omega_a in omega_values:
        omega_a = float(omega_a)
        reference_stability = actuator_stability(reference_rt, omega_a)
        # Two distinct gain studies: fixed absolute g and g/g_crit held fixed while L varies.
        studies = [('fraction', float(fraction)) for fraction in ac['gain_fractions']]
        studies += [('fixed', float(g)) for g in p['comparison_g_per_s']]
        for mode, setting in studies:
            def gain_at(length):
                _, rt, _ = delay_and_boundary(length, c['physics'])
                st = actuator_stability(rt, omega_a)
                return (setting * st['g_crit_per_s'] if mode == 'fraction' else setting), st

            def callback(length):
                g, _ = gain_at(length)
                return evaluate(length, 'classical', g, omega_a)

            rows = []
            for length in lengths:
                length = float(length)
                g, st = gain_at(length)
                if g >= st['g_crit_per_s']:
                    row = dict(arm_length_km=length, scheme='classical', g_per_s=g, omega_a_rad_s=omega_a,
                               operating_status='unstable: no spectrum or key rate evaluated',
                               sns_raw_bps=None, cal_raw_bps=None, sns_bps=None, cal_bps=None)
                else:
                    row = dict(evaluate(length, 'classical', g, omega_a))
                row.update(gain_mode=mode, gain_setting=setting, g_crit_per_s=st['g_crit_per_s'])
                rows.append(row)
                comparisons.append(row)
            reach = roots_for_curve(rows, callback)
            stability_arm = None
            if mode == 'fixed':
                _, first_st = gain_at(float(lengths[0]))
                _, last_st = gain_at(float(lengths[-1]))
                if first_st['g_crit_per_s'] > setting >= last_st['g_crit_per_s']:
                    stability_arm = float(brentq(lambda length: gain_at(length)[1]['g_crit_per_s'] - setting,
                                                 lengths[0], lengths[-1], **root_kwargs()))
            record = dict(omega_a_rad_s=omega_a, gain_mode=mode, gain_setting=setting,
                          reference_g_per_s=(setting * reference_stability['g_crit_per_s'] if mode == 'fraction' else setting),
                          reference_g_crit_per_s=reference_stability['g_crit_per_s'],
                          stability_boundary_arm_km=stability_arm, **reach)
            # A key zero may lie beyond the last coarse stable sample but still before instability.
            # Try the independently source-computed sigma-limit key zero; never set an unstable rate to zero.
            for protocol in ('sns', 'cal'):
                if record[protocol + '_zero_arm_km'] is None:
                    candidate = phase_limit_key_zeros[protocol]
                    if lengths[0] <= candidate <= lengths[-1]:
                        g, st = gain_at(candidate)
                        if g < st['g_crit_per_s']:
                            point = callback(candidate)
                            if point['operating_status'] == 'threshold':
                                record[protocol + '_zero_arm_km'] = candidate
                                record[protocol + '_crossings_arm_km'] = [candidate]
            reaches.append(record)
            if mode == 'fraction':
                g, st = gain_at(reference_length)
                reference_points.append(dict(evaluate(reference_length, 'classical', g, omega_a),
                                             gain_mode=mode, gain_setting=setting, g_crit_per_s=st['g_crit_per_s']))
            print('omega_a [rad/s]', omega_a, mode, setting,
                  'SNS/CAL key-zero arms [km]', record['sns_zero_arm_km'], record['cal_zero_arm_km'],
                  'stability boundary [km]', stability_arm, flush=True)
        # Checkpoint all completed studies; a long scan can be reviewed before it finishes.
        (OUT / 'actuator_checkpoint.json').write_text(json.dumps(dict(config=c, completed_reaches=reaches), indent=2))

    csv_write('actuator_keyrate_curves.csv', comparisons)
    csv_write('actuator_dual_keyrate.csv', dual_rows)
    csv_write('actuator_key_reach.csv', reaches)
    csv_write('actuator_operating_convergence.csv', convergence)
    low_frequency_checks = low_frequency_diagnostics(c, omega_values)
    csv_write('actuator_low_frequency_limit.csv', low_frequency_checks)
    csv_write('actuator_reference_points.csv', reference_points)
    for row in gain_rows:
        row.update(evaluate(reference_length, 'classical', row['g_per_s'], row['omega_a_rad_s']))
    csv_write('actuator_gain_scan.csv', gain_rows)
    # Compare the same PSD against independent spatial quadrature at reference length and safe gain.
    from tfqkd.classical import remote_power_ratio
    quadrature_checks = []
    for omega_a in omega_values:
        st = actuator_stability(reference_rt, float(omega_a))
        g = ac['gain_fractions'][0] * st['g_crit_per_s']
        tau = reference_rt / 2
        f = np.geomspace(1 / c['operation']['tau_max_s'], st['critical_frequency_hz'] * ac['quadrature_frequency_multiple'], ac['quadrature_frequency_points'])
        controller = actuator_controller(g, float(omega_a))
        analytic = remote_controller_ratio(f, tau, controller)
        spatial = remote_power_ratio(f, tau, controller, p['spatial_nodes'])
        quadrature_checks.append(dict(omega_a_rad_s=float(omega_a), max_relative_error=float(np.max(abs(analytic / spatial - 1)))))
    csv_write('actuator_operating_convergence.csv', convergence)
    accepted = [row for row in convergence if max(row['time_relative_change'], row['cap_relative_change'], row['operating_variance_relative_change']) <= p['quadrature_rtol']]
    evidence = dict(config=c, stability_reference=[row for row in stability_rows if row['arm_length_km'] == reference_length],
                    crossing_checks=crossing_checks, reach=reaches, dual_reach=dual_reach,
                    phase_limit_key_zeros=phase_limit_key_zeros, reference_points=reference_points,
                    spatial_quadrature_checks=quadrature_checks,
                    stable_low_frequency_limit_checks=low_frequency_checks,
                    convergence_checks=len(convergence), converged_operating_points=len(cache),
                    max_accepted_time_change=max(row['time_relative_change'] for row in accepted),
                    max_accepted_cap_change=max(row['cap_relative_change'] for row in accepted),
                    max_accepted_variance_change=max(row['operating_variance_relative_change'] for row in accepted),
                    fallback_refinements=sum(row['level'] > 1 for row in convergence))
    (OUT / 'actuator_evidence.json').write_text(json.dumps(evidence, indent=2))
    generate_plots(c, omega_values, stability_rows, comparisons, dual_rows, reaches, dual_reach)
    print('Finished actuator scan:', len(cache), 'converged operating points;', len(reaches), 'reach studies', flush=True)


def low_frequency_diagnostics(config, omega_values):
    p, ac = config['classical'], config['classical_actuator']
    tau, rt, _ = delay_and_boundary(p['length_km'], config['physics'])
    # Williams A9 formal limit; chosen low frequency is the same configured diagnostic as the legacy scan.
    frequency = phase_resonance_frequency(rt) * p['low_frequency_fraction']
    target = float(high_gain_ratio(np.array([frequency]), tau, p['spatial_nodes'])[0])
    rows = []
    for omega_a in omega_values:
        st = actuator_stability(rt, float(omega_a))
        for fraction in ac['gain_fractions']:
            gain = fraction * st['g_crit_per_s']
            # Williams A8 with user pole: compare ONLY finite stable gains, not an inadmissible g->infinity.
            ratio = float(remote_controller_ratio([frequency], tau, actuator_controller(gain, float(omega_a)))[0])
            rows.append(dict(omega_a_rad_s=float(omega_a), g_fraction=fraction, g_per_s=gain,
                             frequency_hz=frequency, remote_ratio=ratio, formal_A9_ratio=target,
                             relative_deviation_from_A9=ratio / target - 1))
    return rows


def generate_plots(config, omega_values, stability_rows, comparisons, dual_rows, reaches, dual_reach=None):
    length = config['classical']['length_km']
    tau, rt, _ = delay_and_boundary(length, config['physics'])
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    spectra = {}
    for index, omega_a in enumerate((omega_values[0], omega_values[-1])):
        st = actuator_stability(rt, float(omega_a))
        for fraction in (config['classical_actuator']['gain_fractions'][1], config['classical_actuator']['gain_fractions'][-2], config['classical_actuator']['gain_fractions'][-1]):
            gain = fraction * st['g_crit_per_s']
            f, _ = actuator_grid(config, length, gain, float(omega_a), level=1)
            f = np.unique(np.r_[f, np.geomspace(config['grid']['f_min_hz'], config['grid']['f_max_hz'], config['grid']['points'])])
            sens = actuator_local_sensitivity(f, rt, gain, float(omega_a))
            remote = remote_controller_ratio(f, tau, actuator_controller(gain, float(omega_a)))
            label = f'{omega_a/(2*np.pi):g} Hz, g/gcrit={fraction:g}'
            axs[0].loglog(f, sens, label=label); axs[1].loglog(f, remote, label=label)
            spectra[f'f_{index}_{fraction}'] = f
            spectra[f'sensitivity_{index}_{fraction}'] = sens
            spectra[f'remote_ratio_{index}_{fraction}'] = remote
        for ax in axs:
            ax.axvline(st['critical_frequency_hz'], color=f'C{index}', linestyle=':',
                       label=f'Marginal frequency, actuator {omega_a/(2*np.pi):g} Hz')
    for ax in axs: ax.axhline(1, color='gray', linewidth=.5); ax.grid(True); ax.legend(fontsize='small')
    axs[0].set_ylabel('Local |1/(1+G)|²'); axs[1].set_ylabel('Remote S_out / S_free (A8)'); axs[1].set_xlabel('Frequency [Hz]')
    fig.tight_layout(); fig.savefig(OUT / 'actuator_suppression.png'); plt.close(fig)
    np.savez_compressed(OUT / 'actuator_spectra.npz', **spectra)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for omega_a in omega_values:
        rows = [r for r in stability_rows if r['omega_a_rad_s'] == omega_a]
        rows.sort(key=lambda r: r['arm_length_km'])
        label = f'Actuator {omega_a/(2*np.pi):g} Hz'
        axs[0].loglog([r['arm_length_km'] for r in rows], [r['g_crit_per_s'] for r in rows], label=label)
        axs[1].loglog([r['arm_length_km'] for r in rows], [r['critical_frequency_hz'] for r in rows], label=label)
    axs[0].set(xlabel='Arm length [km]', ylabel='Critical gain [s⁻¹]')
    axs[1].set(xlabel='Arm length [km]', ylabel='Marginal frequency [Hz]')
    for ax in axs: ax.grid(True); ax.legend(fontsize='small')
    fig.tight_layout(); fig.savefig(OUT / 'actuator_stability.png'); plt.close(fig)
    fractions = config['classical_actuator']['gain_fractions']
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for omega_a in (omega_values[0], omega_values[-1]):
        for fraction in (fractions[0], fractions[-1]):
            rows = [r for r in comparisons if r['omega_a_rad_s'] == omega_a and r['gain_mode'] == 'fraction' and r['gain_setting'] == fraction]
            label = f'{omega_a/(2*np.pi):g} Hz, g/gcrit={fraction:g}'
            x = [r['arm_length_km'] for r in rows]
            axs[0,0].semilogy(x, [r['tau_q_s'] for r in rows], label=label)
            axs[0,1].semilogy(x, [r['duty'] for r in rows], label=label)
            for ax, field in [(axs[1,0], 'sns_bps'), (axs[1,1], 'cal_bps')]:
                ax.semilogy(x, [r[field] if r[field] > 0 else np.nan for r in rows], label=label)
    x = [r['arm_length_km'] for r in dual_rows]
    for ax, field in [(axs[0,0], 'tau_q_s'), (axs[0,1], 'duty'), (axs[1,0], 'sns_bps'), (axs[1,1], 'cal_bps')]:
        ax.semilogy(x, [r[field] if r[field] > 0 else np.nan for r in dual_rows], 'k--', label='Dual-band')
        ax.grid(True); ax.legend(fontsize='x-small'); ax.set_xlabel('Arm length [km]; total distance 2L')
    axs[0,0].set_ylabel('Solved tau_Q [s]'); axs[0,1].set_ylabel('Duty cycle')
    axs[1,0].set_ylabel('SNS-AOPP [bit/s]'); axs[1,1].set_ylabel('CAL [bit/s]')
    if dual_reach is not None:
        for ax, protocol in [(axs[1,0], 'sns'), (axs[1,1], 'cal')]:
            classical_root = next(r[protocol+'_zero_arm_km'] for r in reaches if r['gain_mode']=='fraction')
            ax.axvline(classical_root, color='gray', linestyle=':', label='Classical key zero')
            ax.axvline(dual_reach[protocol+'_zero_arm_km'], color='black', linestyle=':', label='Dual-band key zero')
            ax.legend(fontsize='x-small')
    fig.tight_layout(); fig.savefig(OUT / 'actuator_classical_vs_dual.png'); plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for mode, linestyle in [('fraction', '-'), ('fixed', '--')]:
        settings = fractions if mode == 'fraction' else config['classical']['comparison_g_per_s']
        for setting in settings:
            rows = [r for r in reaches if r['gain_mode'] == mode and r['gain_setting'] == setting]
            x = [r['omega_a_rad_s'] / (2*np.pi) for r in rows]
            for ax, protocol in [(axs[0], 'sns'), (axs[1], 'cal')]:
                y = [r[protocol+'_zero_arm_km'] if r[protocol+'_zero_arm_km'] is not None else np.nan for r in rows]
                ax.semilogx(x, y, linestyle, label=f'{mode} {setting:g}')
    for ax, label in zip(axs, ['SNS-AOPP', 'CAL']):
        ax.set(xlabel='Actuator pole [Hz]', ylabel=f'{label} key-zero arm length [km]');ax.grid(True);ax.legend(fontsize='x-small', ncol=2)
    if dual_reach is not None:
        for ax, protocol in [(axs[0], 'sns'), (axs[1], 'cal')]:
            ax.axhline(dual_reach[protocol+'_zero_arm_km'], color='black', linestyle=':', label='Dual-band')
            ax.legend(fontsize='x-small', ncol=2)
    fig.tight_layout(); fig.savefig(OUT / 'actuator_key_reach.png'); plt.close(fig)


if __name__ == '__main__':
    main()
