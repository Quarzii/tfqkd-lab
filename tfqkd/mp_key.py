"""Zeng mode-pairing asymptotic key per quantum round; independent of TF key code.

Only base IID adjacent-click pairing. No async vacuum-key term or click filter.
Single-photon privacy comes from source decoy bounds, never from removing 25%.
"""
import argparse
import json
from pathlib import Path
import tomllib
import warnings

import numpy as np
from scipy.special import xlogy

from .mp_phase import (IntegrationGrid, components_from_config,
                       phase_difference_variance, zhang_x_error, discrete_average)

COIL_LIMIT = (
    "Published mode-pairing drift parameters were obtained on laboratory fiber spools; "
    "they are not transferable to field routes. For Zhou's four empirical drift values, "
    "alpha=0.5 is rejected within the conditional residual-based 95% log-regression "
    "interval [0.743,1.592], not by a test using published measurement uncertainties. "
    "The data do not identify spatial correlation or a universal spool-to-field coefficient. "
    "Sources: Zhang2025 PDF p5; Zhu2023 Fig2/p3 and Appendix B2/p8; "
    "Zhou async PDF p3-4,26-27; executed fit: results/mp_spatial_analysis/details.json."
)


def _number(name, value, lower=0., upper=None, strict_lower=False):
    value = float(value)
    if not np.isfinite(value) or (value <= lower if strict_lower else value < lower) or (upper is not None and value > upper):
        raise ValueError(f"{name}: finite value {'>' if strict_lower else '>='} {lower} required"
                         + (f", <= {upper}" if upper is not None else ""))
    return value


def _count(name, value):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name}: positive integer required")
    return int(value)


def binary_entropy(probability):
    p = np.asarray(probability, dtype=float)
    if np.any(~np.isfinite(p)) or np.any(p < 0) or np.any(p > 1):
        raise ValueError("Entropy probabilities must lie in [0,1]")
    # Eq. 1 and definition immediately below, zhang2025: binary entropy, including endpoint limits.
    return (-xlogy(p, p) - xlogy(1 - p, 1 - p)) / np.log(2)


def pairing_rate(click_probability, max_gap):
    p = _number("single-click probability", click_probability, upper=1.)
    gap = _count("L_max [round intervals]", max_gap)
    if p == 0:
        return 0.
    if p == 1:
        return .5  # Eq. 4, zeng2022, endpoint with a click in every round.
    # Eq. 4, zeng2022 / Supplement Eq. 72: numerically stable success probability and inverse denominator.
    success = -np.expm1(gap * np.log1p(-p))
    return float(p * success / (1 + success))


def iid_pair_intervals(click_probability, max_gap, active_round_clock_hz):
    p = _number("single-click probability", click_probability, upper=1., strict_lower=True)
    gap = _count("L_max [round intervals]", max_gap)
    clock = _number("active quantum-round clock [Hz]", active_round_clock_hz, strict_lower=True)
    indices = np.arange(1, gap + 1)
    if p == 1:
        weights = np.zeros(gap)
        weights[0] = 1.
    else:
        # Eq. 15/16 and Box 2, zeng2022: geometric waiting times, conditioned on successful gap <= L_max.
        log_weights = (indices - 1) * np.log1p(-p)
        weights = np.exp(log_weights)
        weights /= np.sum(weights)
    # Eq. D5, zhang2025, following paragraph: use discrete pair lengths, not a midpoint quadrature.
    return indices / clock, weights


def phase_average(config, p, max_gap, clock):
    spec = config["phase"]
    distribution = config["intervals"]
    if distribution["mode"] == "iid_zeng":
        times, weights = iid_pair_intervals(p, max_gap, clock)
        interval_note = "Zeng IID geometric gaps; active uniform quantum rounds; reference-frame gaps and L_min>1 absent"
    elif distribution["mode"] == "measured":
        times = np.asarray(distribution["delta_t_s"], dtype=float)
        weights = np.asarray(distribution["counts"], dtype=float)
        if times.ndim != 1 or np.any(~np.isfinite(times)) or np.any(times < 0):
            raise ValueError("Measured pair intervals must be a finite nonnegative vector [s]")
        interval_note = "Explicit measured intervals and counts; no interval density inferred"
    else:
        raise ValueError("Interval mode must be iid_zeng or measured")
    if spec["mode"] == "published_drift":
        drift = _number("sigma_L [rad/s]", spec["sigma_L_rad_s"])
        linewidth = _number("linewidth of each independent laser [Hz]", spec["linewidth_hz"])
        # Eq. 3 / paragraph after Eq. B10, zhang2025: two identical Lorentzian independent lasers.
        variance = 4 * np.pi * linewidth * times + drift**2 * times**2
        phase_metadata = {"model": "Published Gaussian drift + equal Lorentzian lasers; perfect tracking idealization",
                          "sigma_L_rad_s": drift, "linewidth_each_hz": linewidth}
    elif spec["mode"] == "spectral":
        result = phase_difference_variance(components_from_config(config), times,
                                          IntegrationGrid(**config["grid"]), spec.get("T_track_s"))
        variance = result["variance_rad2"]
        phase_metadata = {key: result[key] for key in ("cut_components", "lower_bounds_hz", "T_track_s", "upper_hz")}
    else:
        raise ValueError("Phase mode must be spectral or published_drift; no missing drift is substituted")
    slices = _count("phase slices M", spec["phase_slices"])
    ideal_error = zhang_x_error(variance, slices)
    ed = _number("intrinsic X misalignment", spec["intrinsic_x_error"], upper=.5)
    # Eq. D2, zhang2025: source extension of B10 by the explicitly given intrinsic misalignment factor.
    errors = .5 - (1 - 2 * ed) * (.5 - ideal_error)
    average = float(discrete_average(errors, weights))
    return {"raw_x_error": average, "phase_mode": spec["mode"], "phase_slices": slices,
            "intrinsic_x_error": ed, "distribution": interval_note, "discrete_interval_count": len(times),
            "variance_rad2_min": float(np.min(variance)), "variance_rad2_max": float(np.max(variance)),
            "phase_metadata": phase_metadata}


def _matrix(name, data):
    array = np.asarray(data, dtype=float)
    if array.shape != (3, 3) or np.any(~np.isfinite(array)) or np.any(array < 0):
        raise ValueError(f"{name}: finite nonnegative 3x3 matrix in vacuum/weak/signal order required")
    # Eq. 69/70, zeng2022 supplement: common pair normalizer can rescale these gains/yields; no unit bound invented.
    return array


def decoy_yield_lower(gain, signal, weak):
    q = _matrix("decoy gains", gain)
    m, n = float(signal), float(weak)
    if not np.isfinite(m) or not np.isfinite(n) or not 0 < n < m:
        raise ValueError("Decoy intensities must satisfy signal > weak > 0")
    # Eq. C4/C5, zhang2025: asymptotic limit sets the upper/lower expectation arguments equal to supplied gains.
    weak_term = np.exp(2 * n) * q[1, 1] + q[0, 0] - np.exp(n) * (q[1, 0] + q[0, 1])
    strong_term = np.exp(2 * m) * q[2, 2] + q[0, 0] - np.exp(m) * (q[2, 0] + q[0, 2])
    return float((m**3 * weak_term - n**3 * strong_term) / (m**2 * n**2 * (m - n)))


def decoy_estimates(mu, nu, q_z, q_x, qe_x):
    qz, qx, qe = (_matrix(name, values) for name, values in
                  (("Q_Z", q_z), ("Q_X", q_x), ("QE_X", qe_x)))
    if np.any(qe > qx):
        raise ValueError("Every X error gain must be <= its detection gain")
    yz_raw = decoy_yield_lower(qz, mu, nu)
    # Eq. C5, zhang2025: X-basis total intensities per user over two modes are 2*mu and 2*nu.
    yx_raw = decoy_yield_lower(qx, 2 * mu, 2 * nu)
    # Eq. C6, zhang2025; Eq. 70, zeng2022 supplement: no subtraction of a 25% raw-error floor.
    ye_raw = float((np.exp(4 * nu) * qe[1, 1] + qe[0, 0]
                    - np.exp(2 * nu) * (qe[0, 1] + qe[1, 0])) / (2 * nu)**2)
    if ye_raw < 0:
        raise ValueError(f"Negative decoy error-yield upper bound {ye_raw}: supplied asymptotic gains are inconsistent; no fit or correction applied")
    notes = []
    # Eq. 64, zeng2022 supplement: photon counts are nonnegative; a negative analytical lower bound is vacuous.
    yz, yx = max(0., yz_raw), max(0., yx_raw)
    if qz[2, 2] <= 0:
        raise ValueError("Positive signal Z gain is required to identify its single-photon fraction")
    # Eq. 54/68/69, zeng2022 supplement: Poisson one-photon weights and signal-pair fraction.
    q11 = float(np.exp(-2 * mu) * mu**2 * yz / qz[2, 2])
    if q11 > 1:
        raise ValueError(f"Single-photon fraction lower bound exceeds one ({q11}); check normalization/source data")
    # Eq. C7, zhang2025, asymptotic theta_X=0; Eq. 7, zeng2022, worst-case privacy entropy at error 1/2.
    error_upper_raw = ye_raw / yx if yx > 0 else None
    # Eq. C7, zhang2025: retain the actual error upper bound separately from privacy's worst entropy.
    error_upper = min(1., error_upper_raw) if error_upper_raw is not None else 1.
    e11 = min(.5, error_upper)
    if yz <= 0 or yx <= 0:
        notes.append("VACUOUS_ANALYTIC_DECOY_BOUND: C4/C5 does not establish a positive one-photon yield")
    if yx > 0 and ye_raw / yx >= .5:
        notes.append("NO_SINGLE_PHOTON_PRIVACY_BOUND: upper error interval contains 1/2; privacy contribution is zero")
    return {"q_11_lower": q11, "e_11_x_upper": error_upper,
            "e_11_x_upper_raw": error_upper_raw, "e_11_x_privacy_worst_case": e11, "Y_Z_11_lower_raw": yz_raw,
            "Y_X_11_lower_raw": yx_raw, "YE_X_11_upper": ye_raw, "notes": notes,
            "method": "Zhang C4-C7 asymptotic analytical bounds; Zeng Supplement Note3/Eq70; no finite-key inference"}


def key_rate(click_probability, max_gap, z_pair_fraction, q_11, e_11_x, e_z, f_ec):
    rp = pairing_rate(click_probability, max_gap)
    rs = _number("r_s", z_pair_fraction, upper=1.)
    q = _number("q_(1,1)", q_11, upper=1.)
    ex = _number("single-photon phase-error upper bound", e_11_x, upper=.5)
    ez = _number("e_Z", e_z, upper=.5)
    correction = _number("f_EC", f_ec, lower=1.)
    # Eq. 7, zeng2022: bits per emitted quantum round (two users emit in that round).
    signed = float(rp * rs * (q * (1 - binary_entropy(ex)) - correction * binary_entropy(ez)))
    return {"r_p": rp, "r_s": rs, "signed_bits_per_quantum_round": signed,
            "bits_per_quantum_round": max(0., signed),  # Eq. 7: a nonpositive bound establishes no distillable key.
            "bits_per_potential_pair": 2 * max(0., signed),  # Eq. 1 denominator paragraph, zhu2023: N_pair=N_rounds/2.
            "status": "positive_asymptotic_bound" if signed > 0 else "no_positive_bound"}


def rate_per_second(bits_per_quantum_round, timing):
    clock = _number("active quantum-round clock [Hz]", timing["active_round_clock_hz"], strict_lower=True)
    reference = _number("reference-slot fraction", timing["reference_fraction"], upper=1.)
    recovery = _number("recovery-slot fraction", timing["recovery_fraction"], upper=1.)
    if reference + recovery >= 1:
        raise ValueError("Reference and recovery fractions must leave a positive quantum duty fraction")
    # Eq. 1 normalization and Table III rate columns, zhang2025; Table VI, zhu2023: quantum-round frequency after slots.
    effective = clock * (1 - reference - recovery)
    return {"effective_quantum_rounds_hz": effective, "active_round_clock_hz": clock,
            "reference_fraction": reference, "recovery_fraction": recovery,
            "key_bps": float(bits_per_quantum_round * effective),
            "denominator": "one quantum emission round shared by Alice and Bob; not one successful pair"}


def calculate(config):
    protocol = config["protocol"]
    p, gap = protocol["click_probability"], protocol["max_pairing_gap"]
    phase = phase_average(config, p, gap, config["timing"]["active_round_clock_hz"])
    decoy = config["decoy"]
    qe = _matrix("QE_X", decoy["QE_X"]).copy()
    if decoy["error_mode"] == "phase_weak_equal":
        # Eq. 69, zeng2022 supplement: error gain = gain * raw error; use B10/D2 only for the balanced weak-weak setting.
        qe[1, 1] = _matrix("Q_X", decoy["Q_X"])[1, 1] * phase["raw_x_error"]
        error_note = "Balanced weak X error gain comes from phase model; other error gains remain explicit inputs"
    elif decoy["error_mode"] == "observed":
        error_note = "Privacy bound uses observed error gains; phase calculation is a separate diagnostic"
    else:
        raise ValueError("Decoy error_mode must be phase_weak_equal or observed")
    estimate = decoy_estimates(protocol["mu"], protocol["nu"], decoy["Q_Z"], decoy["Q_X"], qe)
    rate = key_rate(p, gap, protocol["z_pair_fraction"], estimate["q_11_lower"],
                    estimate["e_11_x_privacy_worst_case"], protocol["e_z"], protocol["f_ec"])
    timing = rate_per_second(rate["bits_per_quantum_round"], config["timing"])
    return {"scope": "Zeng2022 base asymptotic MP-QKD, not async Zhou", **rate, **timing,
            "decoy": estimate, "phase": phase, "phase_to_privacy": error_note,
            "limitations": [COIL_LIMIT,
                            "IID Eq4 omits minimum pairing gap, reference-frame gaps, detector dead time and afterpulses.",
                            "B10 is the weak balanced coherent-state ideal-detector phase model; it is not a single-photon error formula.",
                            "Analytical decoy bounds can be looser than the published finite-key LP estimates; asymptotic does not guarantee a larger bound across different estimators.",
                            "No finite-key penalties, phase-tracking estimator model, or async vacuum contribution/click filtering."]
                            + estimate["notes"] + config.get("limitations", [])}


def markdown(result):
    phase, decoy = result.get("phase"), result["decoy"]
    lines = ["# Mode-pairing asymptotic calculation", "", str(result["scope"]), "",
             f"Key rate: {result['key_bps']:.12g} bit/s; {result['bits_per_quantum_round']:.12g} bit/quantum round.",
             f"Potential-pair normalization (Zhu): {result['bits_per_potential_pair']:.12g} bit/(N_rounds/2).", "",
             f"Quantum-round frequency after reference/recovery slots: {result['effective_quantum_rounds_hz']:.12g} Hz.",
             f"r_p = {result['r_p']:.12g}; r_s = {result['r_s']:.12g}.",
             f"Decoy q_(1,1) lower = {decoy['q_11_lower']:.12g}; e_(1,1)^X upper = {decoy['e_11_x_upper']:.12g}.", ""]
    lines += [f"Error used for worst-case privacy entropy = {decoy['e_11_x_privacy_worst_case']:.12g}.", ""]
    if phase:
        lines += [f"Phase mode: {phase['phase_mode']}; discrete intervals: {phase['discrete_interval_count']}.",
                  f"Raw coherent X error = {phase['raw_x_error']:.12g}; {result['phase_to_privacy']}.",
                  phase["distribution"], ""]
    lines += ["## Applicability and limitations", ""] + ["- " + note for note in result["limitations"]]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path, required=True, help="Directory for report.md and result.json")
    args = parser.parse_args()
    with args.config.open("rb") as stream:
        config = tomllib.load(stream)
    result = calculate(config)
    for note in result["decoy"]["notes"]:
        warnings.warn(note)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    (args.output / "report.md").write_text(markdown(result))
    print(json.dumps({"key_bps": result["key_bps"], "status": result["status"], "report": str(args.output / 'report.md')}, indent=2))


if __name__ == "__main__":
    main()
