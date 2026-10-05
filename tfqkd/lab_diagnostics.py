"""Output diagnostics above the unchanged Williams/Bertaina core.

A stability crossing is not a brick-wall bandwidth or an A8 lower bound.
Run: python -m tfqkd.lab_diagnostics INPUT.toml
"""
from __future__ import annotations
import json
import numpy as np
from scipy.optimize import brentq
from .classical import actuator_loop_gain, actuator_stability
from .integration import frequency_grid, PhaseIntegral
from .spectra import free_fiber


def local_suppression_edge(roundtrip_s, gain, omega_a, critical_hz):
    """First unity-sensitivity crossing below the first stability crossing.

    At zero gain sensitivity is identically one: there is no suppressed band.
    At gcrit this evaluates a formal marginal limit, not a runnable controller.
    """
    if gain == 0:
        return 0.0
    def condition(frequency):
        # Eq. A7, williams2008: |1/(1+G)|=1 iff |1+G|^2=1.
        loop = actuator_loop_gain([frequency], roundtrip_s, gain, omega_a)[0]
        return float(abs(1 + loop)**2 - 1)
    return float(brentq(condition, np.nextafter(0., 1.)**0.25 * critical_hz,
                        critical_hz))


def above_band_diagnostic(lab, loop):
    p = lab.config['physics']; omega = lab.settings['actuator']['omega_a_rad_s']
    f = frequency_grid(lab.config); rows = []
    # Eqs. 4,5,7, bertaina2024: threshold variance and arm weights; det noise excluded here.
    threshold = lab.config['operation']['sigma_limit_rad']**2
    weight = p['K'] if lab.settings['scheme']['lasers'] == 'common' else 1
    for length, gc, critical in zip(lab.arms(), loop['arm_g_crit_per_s'], loop['critical_frequency_hz']):
        # Eq. A6, williams2008 + adopted actuator pole: two-way propagation time.
        rt = 2*p['n']*length/p['c_km_s']
        marginal_edge = local_suppression_edge(rt, gc, omega, critical)
        operating_edge = local_suppression_edge(rt, loop['g_per_s'], omega, critical)
        # Eq. 6 and Eq. 4, bertaina2024: diagnostic FREE arm noise above the stability crossing.
        integral = PhaseIntegral(f, free_fiber(f, length, p))
        variance = float(integral.above(critical))
        above_suppression=float(integral.above(marginal_edge))
        rows.append(dict(arm_length_km=length,stability_boundary_hz=critical,
                         marginal_local_suppression_edge_hz=marginal_edge,
                         operating_local_suppression_edge_hz=operating_edge,
                         free_above_stability_variance_rad2=variance,
                         free_above_marginal_suppression_variance_rad2=above_suppression))
    # Eqs. 4,5,7, bertaina2024: combination of independent arm terms with original weights.
    stability_variance = weight * sum(row['free_above_stability_variance_rad2'] for row in rows)
    variance = weight * sum(row['free_above_marginal_suppression_variance_rad2'] for row in rows)
    ratio = variance / threshold
    message = None
    if loop['g_per_s'] == 0 and ratio > 1:
        message = (f'Free fiber noise above the marginal local suppression boundaries contributes {ratio:.6g} '
                   'times the phase-variance threshold. This is a bandwidth diagnostic, '
                   'not a lower bound on the Williams A8 residual. With feedback disabled, '
                   'windows whose lower cutoff includes all these frequencies cannot reach '
                   'the threshold; shorter acquisition windows can. ')
        if loop['mode'].startswith('automatic'):
            message += ('The numerical gain search selected g = 0; no resolved key-rate '
                        f'improvement was found at relative search tolerance '
                        f'{lab.settings["numerics"]["gain_rate_rtol"]:g}. Compare dual-band stabilization.')
        else:
            message += 'Feedback is disabled by the manual gain input.'
    return dict(arms=rows,f_b_hz=min(row['marginal_local_suppression_edge_hz'] for row in rows),
                f_b_definition='first local unity-sensitivity crossing at the formal gcrit limit; not the operating bandwidth or the remote A8 PSD ratio',
                stability_boundary_hz=min(row['stability_boundary_hz'] for row in rows),
                free_above_stability_variance_rad2=stability_variance,
                above_stability_threshold_ratio=stability_variance/threshold,
                operating_suppression_bandwidth_hz=min(row['operating_local_suppression_edge_hz'] for row in rows),
                arm_weight=weight,free_above_band_variance_rad2=variance,
                threshold_variance_rad2=threshold,above_band_threshold_ratio=ratio,message=message,
                interpretation='diagnostic free spectrum integral; actual tau_Q and key benefit must be calculated separately')


def balanced_arm_crossover(lab, omega_a, bracket_km, boundary='suppression'):
    """Root of the diagnostic ratio=1, not a key-rate reach or benefit boundary."""
    p=lab.config['physics'];f=frequency_grid(lab.config)
    # Eq. 6, bertaina2024: free PSD scales linearly with L; compute a unit-length integral once.
    integral=PhaseIntegral(f,free_fiber(f,1.,p))
    # Eqs. 4,5,7, bertaina2024: balanced two-arm diagnostic, same K and threshold as the input.
    weight=p['K'] if lab.settings['scheme']['lasers']=='common' else 1
    threshold=lab.config['operation']['sigma_limit_rad']**2
    def ratio(length):
        # Eq. A6, williams2008 + adopted actuator pole.
        rt=2*p['n']*length/p['c_km_s']
        st=actuator_stability(rt,omega_a);critical=st['critical_frequency_hz']
        if boundary=='suppression':
            cutoff=local_suppression_edge(rt,st['g_crit_per_s'],omega_a,critical)
        elif boundary=='stability':
            cutoff=critical
        else:
            raise ValueError('boundary must be suppression or stability')
        # Eqs. 4,6, bertaina2024, finite upper integration bound retained.
        return float(2*weight*length*integral.above(cutoff)/threshold)
    length=brentq(lambda value:ratio(value)-1,*bracket_km)
    return dict(omega_a_rad_s=float(omega_a),arm_length_km=float(length),physical_total_length_km=float(2*length),
                ratio_at_root=ratio(length),bracket_km=list(bracket_km),
                criterion=f'free above-{boundary}-boundary variance equals threshold; not key reach or gain optimum')


if __name__=='__main__':
    import argparse
    from .lab_inputs import resolve
    from .lab_engine import calculate
    parser=argparse.ArgumentParser();parser.add_argument('config');args=parser.parse_args()
    print(json.dumps(calculate(resolve(args.config)).get('classical_diagnostic'),indent=2))
