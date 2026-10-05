"""Run source comparisons and write all stage-1 numerical evidence and figures."""

# CLI import bootstrap; no calculation settings are changed.
import sys as _sys
from pathlib import Path as _Path
_PROJECT_ROOT = _Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_PROJECT_ROOT))
ROOT = _PROJECT_ROOT

import csv
import copy
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import subprocess

os.environ.setdefault("MPLCONFIGDIR", "/tmp/kvant-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.signal import welch

from tfqkd.config import ROOT, cli_config
from tfqkd.integration import PhaseIntegral, frequency_grid
from tfqkd.noise import generate, discrete_psd_one_sided, power_law_one_sided, stationary_variance, finite_record_expectations
from tfqkd.protocol import duty_cycle, phase_error
from tfqkd.reference import author_functions
from tfqkd.spectra import components, free_fiber

OUTPUT = ROOT / "results"


def write_csv(name, rows):
    with (OUTPUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def percent_difference(value, reference):
    # Numerical relative-error diagnostic for cited source comparisons; not a physical model.
    return float((value / reference - 1) * 100)


def ensemble_summary(values):
    values = np.asarray(values)
    # Empirical standard error across independent Monte Carlo realizations (reporting choice).
    return {"mean": float(values.mean()), "standard_error": float(values.std(ddof=1) / np.sqrt(len(values)))}


def stochastic_checks(c):
    v = c["validation"]
    rng = np.random.default_rng(v["seed"])
    ratios, power_ratios, variances, mean_squares = [], [], [], []
    mean_welch = None
    for _ in range(v["realizations"]):
        x = generate(rng, v["samples"], v["alpha"], v["input_variance"])
        f, measured = welch(x, fs=v["sample_rate_hz"], window=v["welch_window"],
                            nperseg=v["welch_segment"], noverlap=v["welch_overlap"],
                            detrend=False, scaling="density", return_onesided=True)
        positive = f > 0
        target = discrete_psd_one_sided(f[positive], v["sample_rate_hz"], v["alpha"], v["input_variance"])
        power = power_law_one_sided(f[positive], v["sample_rate_hz"], v["alpha"], v["input_variance"])
        band = (f[positive] >= v["comparison_min_hz"]) & (f[positive] <= v["comparison_max_hz"])
        # Eq. 98/99, kasdin1995: dimensionless measured/analytical PSD diagnostics.
        ratios.append(np.mean(measured[positive][band] / target[band]))
        power_ratios.append(np.mean(measured[positive][band] / power[band]))
        # Eq. 111 and Eq. 38, kasdin1995: centered record variance and zero-mean second moment.
        variances.append(np.var(x))
        mean_squares.append(np.mean(x**2))
        mean_welch = measured if mean_welch is None else mean_welch + measured
    # Monte Carlo arithmetic average; comparison targets are Eqs. 98 and 99, kasdin1995.
    mean_welch /= v["realizations"]
    t1 = {"discrete_psd_ratio": ensemble_summary(ratios), "power_law_ratio": ensemble_summary(power_ratios)}
    # Eq. 98, kasdin1995 integrated over both sides; quadrature excludes endpoint evaluation.
    integral, integral_error = quad(lambda f_: 2 * v["input_variance"] / v["sample_rate_hz"] /
                                   (2 * np.sin(np.pi * f_ / v["sample_rate_hz"])) ** v["alpha"],
                                   0, v["sample_rate_hz"] / 2)
    t2 = {"psd_integral": integral, "quadrature_error": integral_error,
          "stationary_variance_eq111": float(stationary_variance(v["alpha"], v["input_variance"])),
          "centered_variance": ensemble_summary(variances), "zero_mean_second_moment": ensemble_summary(mean_squares)}
    for key in ("centered_variance", "zero_mean_second_moment"):
        t2[key]["difference_percent"] = percent_difference(t2[key]["mean"], integral)
        # Empirical Monte Carlo z score; numerical diagnostic, not a source-derived confidence guarantee.
        t2[key]["difference_in_standard_errors"] = (t2[key]["mean"] - integral) / t2[key]["standard_error"]
    # Test thresholds are engineering choices declared in config.toml.
    t1["pass"] = bool(abs(t1["discrete_psd_ratio"]["mean"] - 1) <= v["confidence_sigma"] * t1["discrete_psd_ratio"]["standard_error"])
    t2["stationary_centered_pass"] = bool(abs(t2["centered_variance"]["difference_in_standard_errors"]) <= v["confidence_sigma"])
    t2["zero_mean_variance_pass"] = bool(abs(t2["zero_mean_second_moment"]["difference_in_standard_errors"]) <= v["confidence_sigma"])
    finite = finite_record_expectations(v["samples"], v["alpha"], v["input_variance"])
    t2["finite_record_expectations_eq38"] = finite
    # Eqs. 37-38, kasdin1995: compare unchanged observations with finite-record expectations.
    t2["finite_centered_difference_in_standard_errors"] = (t2["centered_variance"]["mean"] - finite["centered_variance"]) / t2["centered_variance"]["standard_error"]
    t2["finite_centered_pass"] = bool(abs(t2["finite_centered_difference_in_standard_errors"]) <= v["confidence_sigma"])
    # Eq. 98, kasdin1995; lower cutoff 1/T follows the user's revised T2 criterion.
    record_time = v["samples"] / v["sample_rate_hz"]
    band_integral, band_error = quad(lambda f_: 2 * v["input_variance"] / v["sample_rate_hz"] /
                                    (2 * np.sin(np.pi * f_ / v["sample_rate_hz"])) ** v["alpha"],
                                    1 / record_time, v["sample_rate_hz"] / 2)
    # Monte Carlo diagnostic comparing the Eq. 98 integral with the finite-record estimate.
    band_z = (t2["centered_variance"]["mean"] - band_integral) / t2["centered_variance"]["standard_error"]
    t2["revised"] = {"record_time_s": record_time, "lower_hz": 1 / record_time,
                     "band_integral": band_integral, "quadrature_error": band_error,
                     "difference_in_standard_errors": band_z,
                     "band_pass": bool(abs(band_z) <= v["confidence_sigma"]),
                     "finite_expectation_pass": t2["finite_centered_pass"],
                     "pass": bool(abs(band_z) <= v["confidence_sigma"] and t2["finite_centered_pass"])}
    write_csv("T1_psd.csv", [{"f_hz": float(fi), "welch": float(w), "discrete_psd": float(d), "power_law": float(p)}
                              for fi, w, d, p in zip(f[positive], mean_welch[positive], target, power)])
    write_csv("T2_realizations.csv", [{"realization": i, "variance": var, "mean_square": ms}
                                      for i, (var, ms) in enumerate(zip(variances, mean_squares))])
    fig, ax = plt.subplots()
    ax.loglog(f[positive], mean_welch[positive], label="Welch, ensemble mean")
    ax.loglog(f[positive], target, "--", label="Kasdin Eq. 98, one-sided")
    ax.loglog(f[positive], power, ":", label="Kasdin Eq. 99, power law")
    ax.set(xlabel="Frequency [Hz]", ylabel="PSD [test units²/Hz]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT / "T1_psd.png")
    plt.close(fig)
    return t1, t2


def scenario_checks(c, ref):
    f = frequency_grid(c)
    finer = frequency_grid(c, c["grid"]["refined_points"])
    rows, comparisons, convergence, cumulative_rows, cutoff_rows = [], [], [], [], []
    fig, axs = plt.subplots(len(c["scenarios"]), 2, figsize=(11, 18))
    cumulative_fig, cumulative_ax = plt.subplots()
    times = np.asarray(c["validation"]["comparison_times_s"])
    op = c["operation"]
    for i, s in enumerate(c["scenarios"]):
        parts = components(f, s, c)
        integrals = {key: PhaseIntegral(f, value) for key, value in parts.items()}
        integral = integrals["total"]
        tau, status = integral.operating_time(op["sigma_limit_rad"], op["tau_max_s"])
        reference, reference_laser, reference_fiber = ref["calc_spectra"](s["common"], s["cavity"], s["stabilized"])
        expected_psd = reference(f, op["LB_km"], s["delta_L_km"])
        author_tau, _, _ = ref["calc_sigma_tau"](s["common"], s["cavity"], s["stabilized"], s["delta_L_km"], op["LB_km"], op["sigma_limit_rad"])
        author_tau_capped = min(author_tau, op["tau_max_s"])
        # Eq. 4, bertaina2024, author implementation calc_sigma (right rectangles).
        author_t, author_sigma = ref["calc_sigma"](reference, op["LB_km"], s["delta_L_km"], num=c["grid"]["author_curve_points"])
        variance = float(integral.variance(tau))
        # Numerical relative discrepancies, referenced to the unmodified author functions.
        rows.append({"scenario": s["name"], "tau_s": tau, "author_tau_capped_s": float(author_tau_capped),
                     "tau_difference_percent": percent_difference(tau, author_tau_capped),
                     "variance_rad2": variance, "sigma_rad": np.sqrt(variance),
                     "e_phi": float(phase_error(variance)), "duty": duty_cycle(tau, op["tau_ps_s"]), "status": status,
                     "max_psd_relative_error": float(np.max(np.abs(parts["total"] / expected_psd - 1)))})
        refined = PhaseIntegral(finer, components(finer, s, c)["total"])
        refined_tau, _ = refined.operating_time(op["sigma_limit_rad"], op["tau_max_s"])
        convergence.append({"scenario": s["name"],
                            "variance_max_relative_change": float(np.max(np.abs(integral.variance(times) / refined.variance(times) - 1))),
                            "tau_relative_change": float(abs(tau / refined_tau - 1))})
        for cutoff in c["validation"]["extended_cutoffs_hz"]:
            extended_config = copy.deepcopy(c)
            extended_config["grid"]["f_max_hz"] = cutoff
            extended_f = frequency_grid(extended_config, c["grid"]["refined_points"])
            extended = PhaseIntegral(extended_f, components(extended_f, s, c)["total"])
            cutoff_rows.append({"scenario": s["name"], "upper_hz": cutoff,
                                "variance_at_original_tau": float(extended.variance(tau)),
                                "change_from_1MHz_percent": percent_difference(float(extended.variance(tau)), variance)})
        for t in times:
            # Eq. 4, bertaina2024: compare variance; author plots its square root.
            author_variance = float(np.interp(t, author_t, author_sigma**2))
            own = float(integral.variance(t))
            comparisons.append({"scenario": s["name"], "tau_s": t, "variance": own,
                                "author_variance_1000_points": author_variance,
                                "difference_percent": percent_difference(own, author_variance)})
        # Eq. 4, bertaina2024: plot cumulative integral versus its lower frequency limit.
        cumulative_ax.loglog(f[:-1], integral.cumulative[:-1], label=s["name"])
        for fi, vi in zip(f, integral.cumulative):
            cumulative_rows.append({"scenario": s["name"], "lower_hz": fi, "variance_above_rad2": vi})
        for key, style, color in (("laser", ":", "C0"), ("fiber", "--", "C1"), ("detection", "-.", "C2"), ("total", "-", "C3")):
            if np.any(parts[key] > 0):
                axs[i, 0].loglog(f, parts[key], style, color=color, label=key)
                # Eq. 4, bertaina2024: sigma is square root of variance.
                axs[i, 1].loglog(1 / f[:-1], np.sqrt(integrals[key].cumulative[:-1]), style, color=color, label=key)
        axs[i, 0].loglog(f, expected_psd, color="gray", alpha=.5, linewidth=.6, label="author total")
        axs[i, 1].loglog(author_t, author_sigma, color="gray", alpha=.6, linewidth=.7, label="author 1000 points")
        axs[i, 1].axhline(op["sigma_limit_rad"], color="black", linewidth=.5)
        axs[i, 0].set(xlim=(c["plot"]["figure8_f_min_hz"], c["grid"]["f_max_hz"]),
                      ylim=(c["plot"]["figure8_psd_min"], c["plot"]["figure8_psd_max"]), ylabel=f"{s['name']}: PSD [rad²/Hz]")
        axs[i, 1].set(xlim=(1 / c["grid"]["f_max_hz"], op["tau_max_s"]),
                      ylim=(c["plot"]["figure8_sigma_min"], c["plot"]["figure8_sigma_max"]), ylabel="sigma [rad]")
    axs[1, 0].legend(fontsize="small")
    axs[1, 1].legend(fontsize="small")
    axs[-1, 0].set_xlabel("Frequency [Hz]")
    axs[-1, 1].set_xlabel("Integration time [s]")
    fig.tight_layout()
    fig.savefig(OUTPUT / "figure8.png")
    plt.close(fig)
    cumulative_ax.set(xlabel="Lower integration limit [Hz]", ylabel="Variance above cutoff [rad²]")
    cumulative_ax.legend(title="Scenario")
    cumulative_fig.tight_layout()
    cumulative_fig.savefig(OUTPUT / "cumulative_variance.png")
    plt.close(cumulative_fig)
    write_csv("tableI.csv", rows)
    write_csv("figure8_comparison.csv", comparisons)
    write_csv("grid_convergence.csv", convergence)
    write_csv("cumulative_variance.csv", cumulative_rows)
    write_csv("upper_cutoff_sensitivity.csv", cutoff_rows)
    return rows, comparisons, convergence, cutoff_rows


def contour_checks(c, ref):
    v, op = c["validation"], c["operation"]
    mismatches = np.geomspace(v["imbalance_min_km"], v["imbalance_max_km"], v["imbalance_points"])
    f = frequency_grid(c)
    refined_f = frequency_grid(c, c["grid"]["refined_points"])
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    rows = []
    for common, stabilized, panel in ((True, False, 0), (True, True, 1), (False, False, 2), (False, True, 2)):
        for cavity in (False, True):
            own_curve, author_curve = [], []
            for mismatch in mismatches:
                s = {"common": common, "stabilized": stabilized, "cavity": cavity, "delta_L_km": mismatch}
                integral = PhaseIntegral(f, components(f, s, c)["total"])
                own, status = integral.operating_time(op["sigma_limit_rad"], v["contour_tau_max_s"])
                refined_integral = PhaseIntegral(refined_f, components(refined_f, s, c)["total"])
                refined_tau, _ = refined_integral.operating_time(op["sigma_limit_rad"], v["contour_tau_max_s"])
                reference, _, _ = ref["calc_spectra"](common, cavity, stabilized)
                ts, sigmas = ref["calc_sigma"](reference, op["LB_km"], mismatch, num=c["grid"]["author_curve_points"])
                if sigmas[-1] < op["sigma_limit_rad"]:
                    author = np.nan
                else:
                    author = np.interp(op["sigma_limit_rad"], sigmas, ts)
                own = own if status == "threshold" else np.nan
                own_curve.append(own)
                author_curve.append(author)
                rows.append({"common": common, "cavity": cavity, "stabilized": stabilized, "delta_L_km": mismatch,
                             "tau_s": own, "author_tau_s": author,
                             "difference_percent": percent_difference(own, author),
                             "refined_tau_s": refined_tau,
                             "grid_change_percent": percent_difference(own, refined_tau)})
            label = ("stable laser" if cavity else "free laser") + (", fiber stable" if stabilized else ", fiber free")
            line, = axs[panel].loglog(own_curve, mismatches, label=label)
            axs[panel].loglog(author_curve, mismatches, ":", color=line.get_color(), alpha=.7)
    for ax, title in zip(axs, ("Common laser, free fiber", "Common laser, stable fiber", "Independent lasers")):
        ax.set(xlabel="Threshold time [s]", title=title)
        ax.grid(True)
        ax.legend(fontsize="x-small")
    axs[0].set_ylabel("Length mismatch [km]")
    fig.suptitle("Figure 2: solid = present quadrature; dotted = author 1000-point curves")
    fig.tight_layout()
    fig.savefig(OUTPUT / "figure2.png")
    plt.close(fig)
    write_csv("figure2_comparison.csv", rows)
    worst = max(rows, key=lambda r: abs(r["difference_percent"]) if np.isfinite(r["difference_percent"]) else -np.inf)
    worst_reference, _, _ = ref["calc_spectra"](worst["common"], worst["cavity"], worst["stabilized"])
    author_refinement = []
    integral_refinement = []
    for count in (c["grid"]["author_curve_points"], c["grid"]["author_points"], c["grid"]["points"], c["grid"]["refined_points"]):
        ts, sigmas = ref["calc_sigma"](worst_reference, op["LB_km"], worst["delta_L_km"], num=count)
        author_refined_tau = float(np.interp(op["sigma_limit_rad"], sigmas, ts))
        author_refinement.append({"points": count, "author_tau_s": author_refined_tau,
                                  "difference_percent": percent_difference(worst["refined_tau_s"], author_refined_tau)})
        # Eq. 4, bertaina2024: fixed tau isolates quadrature error from threshold finding.
        local_f = frequency_grid(c, count)
        exact_cutoff = PhaseIntegral(local_f, components(local_f, worst, c)["total"])
        integral_refinement.append({"points": count,
                                   "fixed_tau_s": worst["refined_tau_s"],
                                   "trapezoid_variance": float(exact_cutoff.variance(worst["refined_tau_s"])),
                                   "author_variance": float(np.interp(worst["refined_tau_s"], ts, sigmas**2))})
    write_csv("figure2_integral_convergence.csv", integral_refinement)
    fig, ax = plt.subplots()
    for key, label in (("trapezoid_variance", "Trapezoids, exact lower cutoff"), ("author_variance", "Authors' right rectangles")):
        ax.semilogx([r["points"] for r in integral_refinement], [r[key] for r in integral_refinement], "o-", label=label)
    # Eq. 4, bertaina2024: fixed reference time was solved for this variance threshold.
    ax.axhline(op["sigma_limit_rad"]**2, color="gray", linestyle=":", label="Threshold variance")
    ax.set(xlabel="Logarithmic frequency-grid points", ylabel="Variance at fixed tau [rad²]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT / "figure2_integral_convergence.png")
    plt.close(fig)
    return {"points": len(rows), "max_abs_tau_difference_percent": float(np.nanmax(np.abs([r["difference_percent"] for r in rows]))),
            "max_abs_grid_change_percent": float(np.nanmax(np.abs([r["grid_change_percent"] for r in rows]))),
            "worst_author_difference": worst, "author_refinement_at_worst_point": author_refinement,
            "integral_refinement_at_worst_point": integral_refinement,
            "censored_points": sum(not np.isfinite(r["tau_s"]) for r in rows)}


def measurement_check(c):
    v = c["validation"]
    source = ROOT / "sources/data/free_fiber_meas.txt"
    f, measured = np.loadtxt(source, unpack=True)
    # Eq. 6 and Appendix G, bertaina2024; author notebook Figure 6: measured double pass.
    modeled = v["measurement_roundtrip_factor"] * free_fiber(f, v["measurement_length_km"], c["physics"])
    edges = v["measurement_edges_hz"]
    ranges = list(zip(edges[:-1], edges[1:])) + [(edges[0], edges[-1]), (edges[1], edges[3]),
                                              (float(f.min()), c["physics"]["fc1_hz"])]
    rows = []
    for low, high in ranges:
        mask = (f >= low) & (f < high) & (measured > 0)
        # Eq. 6, bertaina2024: descriptive log-residuals and least-squares slopes only;
        # neither model coefficients nor frequency windows are fitted to the measurements.
        ratio = measured[mask] / modeled[mask]
        x = np.log10(f[mask])
        measured_slope = np.polyfit(x, np.log10(measured[mask]), 1)[0]
        model_slope = np.polyfit(x, np.log10(modeled[mask]), 1)[0]
        # Eq. 6, bertaina2024: low-f asymptote, shown separately without replacing Eq. 6.
        asymptote = v["measurement_roundtrip_factor"] * c["physics"]["l"] * v["measurement_length_km"] / f[mask]**2
        rows.append({"low_hz": low, "high_hz_exclusive": high, "points": int(mask.sum()),
                     "median_measurement_over_model": float(np.median(ratio)),
                     "median_residual_db": float(10 * np.log10(np.median(ratio))),
                     "rms_log10_residual": float(np.sqrt(np.mean(np.log10(ratio)**2))),
                     "measured_slope": float(measured_slope), "model_slope": float(model_slope),
                     "slope_difference": float(measured_slope - model_slope)})
        rows[-1]["strict_region_below_fc1"] = bool(high <= c["physics"]["fc1_hz"])
        rows[-1]["median_measurement_over_low_f_asymptote"] = float(np.median(measured[mask] / asymptote))
        rows[-1]["low_f_asymptote_slope"] = float(np.polyfit(x, np.log10(asymptote), 1)[0])
    write_csv("T4_residuals.csv", rows)
    write_csv("T4_spectrum.csv", [{"f_hz": fi, "measured_rad2_hz": mi, "model_rad2_hz": si}
                                 for fi, mi, si in zip(f, measured, modeled)])
    fig, axs = plt.subplots(2, 1, sharex=True)
    axs[0].loglog(f, measured, alpha=.6, label="Authors' measured 114 km trace")
    axs[0].loglog(f, modeled, label="4 × Eq. 6, fixed coefficients")
    axs[0].set_ylabel("PSD [rad²/Hz]")
    axs[0].legend(fontsize="small")
    # Dimensionless residual diagnostic of Eq. 6, bertaina2024.
    axs[1].loglog(f, measured / modeled)
    axs[1].axhline(1, color="black", linewidth=.5)
    axs[1].set(xlabel="Frequency [Hz]", ylabel="Measured / model", xlim=(edges[0], edges[-1]))
    fig.tight_layout()
    fig.savefig(OUTPUT / "T4_fiber114km.png")
    plt.close(fig)
    fig, ax = plt.subplots()
    below = f < c["physics"]["fc1_hz"]
    # Eq. 6, bertaina2024, double-pass low-frequency asymptote; distinct from full Eq. 6.
    asymptote_all = v["measurement_roundtrip_factor"] * c["physics"]["l"] * v["measurement_length_km"] / f**2
    ax.loglog(f[below], measured[below], alpha=.6, label="Measured")
    ax.loglog(f[below], modeled[below], label="4 × full Eq. 6")
    ax.loglog(f[below], asymptote_all[below], "--", label="4 l L / f² asymptote")
    ax.set(xlabel="Frequency [Hz]", ylabel="PSD [rad²/Hz]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT / "T4_below_cutoff.png")
    plt.close(fig)
    return rows


def main():
    global OUTPUT
    c = cli_config()
    OUTPUT = Path(c.get("_output_dir", OUTPUT))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    t1, t2 = stochastic_checks(c)
    print("T1", json.dumps(t1), flush=True)
    print("T2", json.dumps(t2), flush=True)
    ref = author_functions()
    rows, comparisons, convergence, cutoffs = scenario_checks(c, ref)
    contour = contour_checks(c, ref)
    measurement = measurement_check(c)
    sources = ["sources/data/QKD.ipynb", "sources/papers/bertaina2024.pdf", "sources/papers/kasdin1995.pdf", "sources/papers/clivati2022-1.pdf", "sources/papers/clivati2022-2.pdf", "sources/data/free_fiber_meas.txt"]
    evidence = {"config": c, "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
                "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sources},
                "T1": t1, "T2": t2, "tableI": rows, "figure8": comparisons, "convergence": convergence,
                "figure2": contour, "upper_cutoff_sensitivity": cutoffs, "T4": measurement, "figure3": "deferred_to_stage2_by_user"}
    evidence["acceptance"] = {
        "T1": "pass" if t1["pass"] else "failed",
        "T2": "pass" if t2["revised"]["pass"] else "failed",
        "T2_original_full_integral_diagnostic": "pass" if t2["stationary_centered_pass"] else "failed",
        "T2_finite_record_diagnostic": "pass" if t2["finite_centered_pass"] else "failed",
        "T3": "discretization_result_converged_Figure3_in_stage2",
        "T4": "strict_region_f_below_fc1_residuals_reported_no_uncertainties",
        "stage": "accepted_by_user_with_reported_discrepancies",
    }
    checked = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", str(ROOT / "tests"), "-p", "test_stage1.py", "-v"],
                             cwd=ROOT, capture_output=True, text=True)
    (OUTPUT / "unit_tests.txt").write_text(checked.stdout + checked.stderr)
    evidence["unit_tests_exit_code"] = checked.returncode
    (OUTPUT / "evidence.json").write_text(json.dumps(evidence, indent=2, ensure_ascii=False))
    print("Figure 2:", contour)
    print("Maximum grid change:", max(r["variance_max_relative_change"] for r in convergence))
    print("Evidence written to", OUTPUT)
    print("Stage 1 acceptance:", evidence["acceptance"])
    return 0 if t1["pass"] and t2["revised"]["pass"] and checked.returncode == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
