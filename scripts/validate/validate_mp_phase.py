"""Analytic gate first, then unfitted MP phase/graphic comparisons; no key formula."""
from dataclasses import replace
from pathlib import Path
import argparse
import hashlib
import json
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import sici

from tfqkd.mp_phase import (IntegrationGrid, fiber_component, laser_components,
                           component_variance, phase_difference_variance,
                           fiber_short_interval, equivalent_fiber_drift,
                           zhou_x_error, zhang_x_error, discrete_average)

ROOT = Path(__file__).resolve().parents[2]


def independent_fiber_ratio(interval, cutoff, numerics):
    """Independent adaptive quadrature of the dimensionless Eq. 6 increment."""
    # Eq. 6, bertaina2024 + Eq. 1, didomenico2010: a=fc*dt and R=integral sinc(a*x)^2/(1+x)^2 dx.
    a = cutoff * interval
    epsilon, limit = numerics["quadrature_epsabs"], numerics["quadrature_limit"]
    low = quad(lambda x: np.sinc(a * x)**2 / (1 + x)**2,
               0, 1, epsabs=epsilon, limit=limit)[0]
    # Same Eqs., change x=exp(u): smooth middle frequency interval, no spectrum/grid import.
    middle = quad(lambda u: np.exp(u) * np.sinc(a * np.exp(u))**2 / (1 + np.exp(u))**2,
                  0, np.log(1 / a), epsabs=epsilon, limit=limit)[0]
    # Same Eqs., change z=a*x: convergent oscillatory tail, independently integrated by scipy.quad.
    tail = quad(lambda z: a * np.sinc(z)**2 / (a + z)**2,
                1, numerics["reference_phase_extent"], epsabs=epsilon, limit=limit)[0]
    # Same Eqs., sin^2 <= 1: bound on the omitted reference tail, not a fabricated uncertainty.
    tail_bound = a / (3 * np.pi**2 * numerics["reference_phase_extent"]**3)
    return low + middle + tail, tail_bound


def analytic_checks(config, grid):
    numerics, fiber = config["numerical_validation"], config["fiber_model"]
    length = config["zhou"]["lengths_km"][0]
    term = fiber_component("fiber", dict(fiber, length_km=length))
    times = np.asarray(numerics["analytic_intervals_s"])
    leading = fiber_short_interval(times, fiber["l"], length, fiber["fc1_hz"])
    corrected = fiber_short_interval(times, fiber["l"], length, fiber["fc1_hz"], next_order=True)
    integral, _ = component_variance(term, times, grid)
    reference = [independent_fiber_ratio(t, fiber["fc1_hz"], numerics) for t in times]
    ratios = np.asarray([r[0] for r in reference])
    # Eq. 6, bertaina2024, derived dimensionless ratio: restore physical units to reference integral.
    exact_reference = leading * ratios
    errors = np.abs(integral / exact_reference - 1)  # Eq. 6-derived reference discrepancy, arithmetic.
    if np.max(errors) > numerics["relative_tolerance"]:
        raise AssertionError(f"Independent fiber integral mismatch {np.max(errors)}")
    rows = []
    for index, t in enumerate(times):
        rows.append({"interval_s": t, "fc_times_interval": fiber["fc1_hz"] * t,
                     "grid_rad2": integral[index], "independent_reference_rad2": exact_reference[index],
                     "relative_grid_error": errors[index], "reference_tail_bound_in_ratio": reference[index][1],
                     "leading_rad2": leading[index], "cubic_corrected_rad2": corrected[index],
                     "leading_relative_overestimate": leading[index] / exact_reference[index] - 1,
                     "corrected_relative_error": corrected[index] / exact_reference[index] - 1,
                     "observed_first_correction_coefficient": (1 - ratios[index]) / (fiber["fc1_hz"] * t)})
    # Eq. 6, bertaina2024, derived next-order coefficient; numerical gate from configured accuracy.
    cubic_limit = np.pi**2 / 3
    if abs(rows[0]["observed_first_correction_coefficient"] / cubic_limit - 1) > numerics["curve_relative_grid_tolerance"]:
        raise AssertionError("Small-interval cubic coefficient did not approach the derived value")
    linewidth = numerics["white_linewidth_hz"]
    white = laser_components("laser", {"model": "lorentzian_white_frequency", "linewidth_hz": linewidth})[0]
    white_t = np.asarray(numerics["white_check_intervals_s"])
    white_d, _ = component_variance(white, white_t, grid)
    # Eq. 4 and white-noise paragraph, didomenico2010: closed one-sided finite-band result.
    si, _ = sici(2 * np.pi * grid.upper_hz * white_t)
    finite = 4 * linewidth / np.pi * (np.pi * white_t * si
              - np.sin(np.pi * grid.upper_hz * white_t)**2 / grid.upper_hz)
    # Eq. 1, didomenico2010 + Eq. 3, zhang2025: infinite-band one-laser limit.
    infinite = 2 * np.pi * linewidth * white_t
    white_errors = np.abs(white_d / finite - 1)  # Eq. 4 comparison, arithmetic.
    if np.max(white_errors) > numerics["relative_tolerance"]:
        raise AssertionError(f"White frequency finite-band mismatch {np.max(white_errors)}")
    convergence = []
    for points in (numerics["fast_points"], grid.points, numerics["refined_points"]):
        test_grid = replace(grid, points=points)
        fd, _ = component_variance(term, times, test_grid)
        wd, _ = component_variance(white, white_t, test_grid)
        convergence.append({"points": points,
                            "fiber_max_relative_error_vs_quad": float(np.max(np.abs(fd / exact_reference - 1))),
                            "white_max_relative_error_vs_closed_finite_band": float(np.max(np.abs(wd / finite - 1)))})
    return {"status": "PASS", "performed_before_experimental_comparison": True,
            "relative_tolerance": numerics["relative_tolerance"],
            "fiber_rows": rows, "first_cubic_coefficient_pi2_over3": cubic_limit,
            "white": {"times_s": white_t.tolist(), "grid_rad2": white_d.tolist(),
                      "closed_finite_rad2": finite.tolist(), "infinite_rad2": infinite.tolist(),
                      "max_relative_error": float(np.max(white_errors)),
                      "relative_upper_band_deficit": (1 - finite / infinite).tolist()},
            "grid_convergence": convergence}


def drift_comparison(config):
    fiber, zhou, zhang = config["fiber_model"], config["zhou"], config["zhang"]
    rows = []
    for length, published in zip(zhou["lengths_km"], zhou["sigma_published_rad_s"]):
        predicted = equivalent_fiber_drift(fiber["l"], length, fiber["fc1_hz"])
        # Eq. 6, bertaina2024 + published Fig. 3, zhou2023_async: unfitted order-of-magnitude comparison.
        rows.append({"source": "zhou2023_async", "total_length_km": length,
                     "published_rad_s": published, "predicted_rad_s": predicted,
                     "ratio_predicted_published": predicted / published,
                     "relative_difference_percent": 100 * (predicted / published - 1),
                     "published_over_sqrt_length": published / np.sqrt(length)})
    predicted = equivalent_fiber_drift(fiber["l"], zhang["length_km"], fiber["fc1_hz"])
    rows.append({"source": "zhang2025", "total_length_km": zhang["length_km"],
                 "published_rad_s": zhang["sigma_published_rad_s"], "predicted_rad_s": predicted,
                 "ratio_predicted_published": predicted / zhang["sigma_published_rad_s"],
                 "relative_difference_percent": 100 * (predicted / zhang["sigma_published_rad_s"] - 1)})
    # Eq. 6 derived sqrt(L) law vs published endpoint ratio; this is a diagnostic, not a parameter fit.
    observed_endpoint_exponent = np.log(zhou["sigma_published_rad_s"][-1] / zhou["sigma_published_rad_s"][0]) / np.log(zhou["lengths_km"][-1] / zhou["lengths_km"][0])
    return {"rows": rows, "predicted_length_exponent": .5,
            "published_endpoint_exponent": observed_endpoint_exponent,
            "noise_coefficients_fitted": False,
            "interpretation": "Different line coefficients; order-of-magnitude comparison, not apparatus calibration"}


def spectra_for_experiment(config, family, length, laser=True):
    fiber = config["fiber_model"]
    setup = config[family]
    if family == "zhou":
        index = setup["lengths_km"].index(length)
        arms = setup["arm_lengths_km"][index]
    else:
        arms = setup["arm_lengths_km"]
    # Eq. 6, bertaina2024; Fig. 2, zeng2022: each arm once, independent noise, no K.
    components = [fiber_component(name, dict(fiber, length_km=arm))
                  for name, arm in zip(("fiber_A", "fiber_B"), arms)]
    if laser:
        # Eq. 3, zhang2025 + white-noise paragraph, didomenico2010: declared Lorentzian approximation.
        components += [c for name in ("laser_A", "laser_B") for c in laser_components(name,
                       {"model": "lorentzian_white_frequency",
                        "linewidth_hz": setup["lorentzian_white_frequency_linewidth_per_laser_hz"]})]
    return components


def residual_statistics(predicted, observed):
    # Error-curve comparison only, not a statistical significance test or fitted physical model.
    residual = predicted - observed
    return {"rms_probability": float(np.sqrt(np.mean(residual**2))),
            "mean_probability": float(np.mean(residual)),
            "maximum_absolute_probability": float(np.max(np.abs(residual))),
            "statistical_significance": "NOT_ESTABLISHED; raw covariance/confidence metadata absent"}


def compare_curves(config, grid, digitization, output):
    rows, blocked = [], []
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    colors = ("tab:green", "black", "tab:blue", "tab:red")
    for curve in digitization["curves"]:
        observed = np.loadtxt(ROOT / curve["csv"], delimiter=",", skiprows=1)
        if curve["id"].startswith("zhu"):
            blocked.append({"curve": curve["id"], "csv": curve["csv"],
                            "points": len(observed), "gap_bin": curve["pulse_gap_bin"],
                            "status": "NO_PREDICTION", "missing_inputs": config["zhu"]["missing_inputs"]})
            continue
        times, data = observed[:, 0], observed[:, 1]
        family = "zhou" if curve["id"].startswith("zhou") else "zhang"
        components = spectra_for_experiment(config, family, curve["length_km"])
        result = phase_difference_variance(components, times, grid)
        finer = phase_difference_variance(components, times, replace(grid, points=config["numerical_validation"]["refined_points"]))
        # Eq. 1, didomenico2010: quadrature refinement, not a fitted scale factor.
        relative_grid_change = float(np.max(np.abs(result["variance_rad2"] / finer["variance_rad2"] - 1)))
        if relative_grid_change > config["numerical_validation"]["curve_relative_grid_tolerance"]:
            raise AssertionError(f"Curve phase quadrature not converged: {curve['id']}, {relative_grid_change}")
        if family == "zhou":
            # Eq. 2, zhou2023_async: source reference curve uses its published Gaussian drift rates.
            source_model = zhou_x_error(curve["published_sigma_rad_s"]**2 * times**2,
                                       times, config["zhou"]["visibility_v2"], curve["residual_frequency_hz"])
            predicted = zhou_x_error(result["variance_rad2"], times,
                                    config["zhou"]["visibility_v2"], curve["residual_frequency_hz"])
            axis = axes[0, 0] if curve["id"].startswith("zhou_3a") else axes[0, 1]
            index = int(curve["id"].rsplit("_", 1)[1])
            label = f"{curve['residual_frequency_hz']:g} Hz" if curve["id"].startswith("zhou_3a") else f"{curve['length_km']:g} km"
            axis.plot(times * 1e6, data, ".", color=colors[index], markersize=3, label=label + " published")
            axis.plot(times * 1e6, predicted, "-", color=colors[index], linewidth=1)
        else:
            # Eq. 3/B10, zhang2025: the paper's stated simplified phase model, with perfect tracking.
            source_d = 4 * np.pi * config["zhang"]["lorentzian_white_frequency_linewidth_per_laser_hz"] * times + curve["published_sigma_rad_s"]**2 * times**2
            source_model = zhang_x_error(source_d, config["zhang"]["phase_slices"])
            predicted = zhang_x_error(result["variance_rad2"], config["zhang"]["phase_slices"])
            axes[1, 0].semilogx(times * 1e6, data, ".", color="tab:blue", markersize=3)
            axes[1, 0].semilogx(times * 1e6, predicted, "-", color="tab:red", linewidth=1)
            axes[1, 0].semilogx(times * 1e6, source_model, "--", color="black", linewidth=1)
        csv = output / (curve["id"] + "_comparison.csv")
        np.savetxt(csv, np.c_[times, data, result["variance_rad2"], source_model, predicted],
                   delimiter=",", header="delta_t_s,published_x_error,transferred_psd_variance_rad2,source_model_x_error,transferred_psd_x_error", comments="")
        rows.append({"curve": curve["id"], "points": len(data), "csv": str(csv.relative_to(ROOT)),
                     "total_length_km": curve["length_km"],
                     "paper_model_against_digitized": residual_statistics(source_model, data),
                     "unfitted_transferred_psd_against_digitized": residual_statistics(predicted, data),
                     "max_relative_grid_change": relative_grid_change,
                     "cut_components": result["cut_components"],
                     "laser_model_assumption": config[family]["laser_model_assumption"],
                     "missing_inputs": config[family]["missing_inputs"],
                     "validation_type": "Conditional Gaussian/Lorentzian source-model comparison; not apparatus-specific PSD validation"})
    axes[0, 0].set_title("Zhou Fig. 3(a): beat frequency, 201.86 km")
    axes[0, 1].set_title("Zhou Fig. 3(b): length; lines use transferred PSD")
    axes[1, 0].set_title("Zhang Fig. 1(b): blue data, black source, red transferred PSD")
    axes[1, 1].axis("off")
    axes[1, 1].text(0, .95, "No fit: l=44 rad² Hz/km, fc1=100 Hz\nfrom another line (Bertaina Table III).\n\nZhu bin averages are digitized separately.\nPrediction blocked: estimator residual and\nwithin-bin interval distribution unavailable.\n\nRaw coherent X errors only; no key formula.",
                    va="top", fontsize=11)
    for axis in axes.ravel()[:3]:
        axis.set_xlabel("Pair interval [µs]")
        axis.set_ylabel("Raw X error probability")
        axis.grid(alpha=.2)
    axes[0, 0].legend(fontsize=7)
    axes[0, 1].legend(fontsize=7)
    figure.savefig(output / "x_error_comparison.png", dpi=180)
    plt.close(figure)
    figure, axis = plt.subplots(figsize=(7, 4), constrained_layout=True)
    for curve in digitization["curves"]:
        if curve["id"].startswith("zhu"):
            data = np.loadtxt(ROOT / curve["csv"], delimiter=",", skiprows=1)
            axis.plot(data[:, 0], data[:, 1], "o-", label=str(curve["pulse_gap_bin"]))
    axis.set_xlabel("Published fiber length [km]")
    axis.set_ylabel("Strong-reference bin-averaged raw X error")
    axis.set_title("Zhu Fig. 6(a): source data only; no unreported tracking input")
    axis.legend(title="Pulse-gap bin", fontsize=8)
    figure.savefig(output / "zhu_digitized.png", dpi=180)
    plt.close(figure)
    return rows, blocked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "sources/data/mp_phase/comparison_inputs.json")
    parser.add_argument("--digitization", type=Path, default=ROOT / "sources/data/mp_phase/digitized/digitization.json")
    parser.add_argument("--output", type=Path, default=ROOT / "results/mp_phase")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    config = json.loads(args.config.read_text())
    grid = IntegrationGrid(**config["grid"])
    analytic = analytic_checks(config, grid)
    # Persist the passing analytic gate before any experimental curve calculation.
    (args.output / "analytic.json").write_text(json.dumps(analytic, indent=2) + "\n")
    print("Analytic gate PASS", flush=True)
    drift = drift_comparison(config)
    digitization = json.loads(args.digitization.read_text())
    curves, blocked = compare_curves(config, grid, digitization, args.output)
    table = config["zhu"]["table_iv_202km"]
    averages = {}
    for kind in ("strong", "qkd"):
        counts, errors = np.asarray(table[kind + "_pairs"]), np.asarray(table[kind + "_error_pairs"])
        # Table IV, zhu2023 + Eq. D5, zhang2025: discrete sum over actual observed counts.
        averages[kind] = {"raw_x_error": discrete_average(errors / counts, counts),
                          "pair_count": int(np.sum(counts)), "error_pair_count": int(np.sum(errors)),
                          "source": "Table IV, zhu2023; no synthetic bin midpoint or within-bin distribution"}
    baseline = args.output / "before_hashes.json"
    if baseline.exists():
        hashes = json.loads(baseline.read_text())
        changed = [path for path, old in hashes.items() if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != old]
        if changed:
            raise AssertionError("Pre-existing TF files changed: " + str(changed))
        unchanged = {"file_count": len(hashes), "changed": changed}
    else:
        unchanged = {"status": "No baseline supplied"}
    result = {"scope": "Mode-pairing phase only; no key implementation", "analytic": analytic,
              "drift_comparison": drift, "curve_comparisons": curves, "blocked_curves": blocked,
              "discrete_observed_averages": averages, "tf_unchanged": unchanged,
              "runtime_s": time.perf_counter() - started,
              "source_inputs": str(args.config.relative_to(ROOT)),
              "digitization": str(args.digitization.relative_to(ROOT)),
              "no_parameter_fitting": True}
    (args.output / "details.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"analytic": "PASS", "compared_curves": len(curves), "blocked_zhu_curves": len(blocked),
                      "drift": drift, "runtime_s": result["runtime_s"]}, indent=2))


if __name__ == "__main__":
    main()
