"""T6 and the source-supported phase-model parts of T7; no invented loop model."""

# CLI import bootstrap; no calculation settings are changed.
import sys as _sys
from pathlib import Path as _Path
_PROJECT_ROOT = _Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_PROJECT_ROOT))
ROOT = _PROJECT_ROOT

import copy
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/kvant-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tfqkd.config import ROOT, cli_config
from tfqkd.integration import frequency_grid, PhaseIntegral
from tfqkd.linewidth import white_linewidth, numerical_linewidth
from tfqkd.spectra import components
from tfqkd.transfers import common_laser_power


def main():
    c = cli_config()
    out = Path(c.get("_output_dir", ROOT / "results"))
    out.mkdir(parents=True, exist_ok=True)
    white = []
    for h0 in c["linewidth"]["white_levels_hz2_per_hz"]:
        analytic = white_linewidth(h0)
        numerical = numerical_linewidth(h0, c["linewidth"]["quadrature_tolerance"])
        # Eq. 5, didomenico2010: numerical relative residual against white-noise limit.
        relative = numerical / analytic - 1
        white.append({"h0_hz2_per_hz": h0, "analytic_hz": analytic, "numerical_hz": numerical,
                      "relative_difference": relative, "pass": abs(relative) <= c["linewidth"]["relative_tolerance"]})
    filtered = []
    for fc in c["linewidth"]["cutoffs_hz"]:
        filtered.append({"fc_hz": fc, "h0_hz2_per_hz": c["linewidth"]["filtered_h0_hz2_per_hz"],
                         "numerical_fwhm_hz": numerical_linewidth(c["linewidth"]["filtered_h0_hz2_per_hz"],
                                                                   c["linewidth"]["quadrature_tolerance"], fc)})
    f = frequency_grid(c)
    n_values = np.linspace(c["sensitivity"]["n_min"], c["sensitivity"]["n_max"], c["sensitivity"]["points"])
    rows = []
    for s in c["scenarios"]:
        baseline_integral = PhaseIntegral(f, components(f, s, c)["total"])
        tau, _ = baseline_integral.operating_time(c["operation"]["sigma_limit_rad"], c["operation"]["tau_max_s"])
        baseline = float(baseline_integral.variance(tau))
        for n in n_values:
            alternate = copy.deepcopy(c)
            alternate["physics"]["n"] = float(n)
            variance = float(PhaseIntegral(f, components(f, s, alternate)["total"]).variance(tau))
            # Eqs. 4, 5, 7, bertaina2024: sensitivity at fixed time, relative to configured n.
            rows.append({"scenario": s["name"], "n": float(n), "fixed_tau_s": tau,
                         "reference_n": c["physics"]["n"], "variance_rad2": variance,
                         "change_percent": (variance / baseline - 1) * 100})
    # Eq. 5, bertaina2024: exact equal-arm limit is zero at every positive frequency.
    equal_arm_max = float(np.max(np.abs(common_laser_power(f, 0, c["physics"]))))
    with (out / "T7_n_sensitivity.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, ax = plt.subplots()
    for s in c["scenarios"]:
        subset = [r for r in rows if r["scenario"] == s["name"]]
        ax.plot([r["n"] for r in subset], [r["change_percent"] for r in subset], label=s["name"])
    ax.set(xlabel="Refractive index n", ylabel="Variance change at fixed tau [%]")
    ax.legend(title="Scenario")
    fig.tight_layout()
    fig.savefig(out / "T7_n_sensitivity.png")
    plt.close(fig)
    evidence = {"config": c, "T6_white": white, "T6_filtered_Fig1": filtered,
                "T7_equal_arm_max": equal_arm_max, "T7_equal_arm_pass": equal_arm_max == 0,
                "T7_n_sensitivity": rows, "T7_loop_limits": "actuator_evidence.json: finite actuator-pole g_crit; g>=g_crit rejected. Ideal C=g/s limits are archived in results/stage2_fixed_time/classical_evidence.json and are not apparatus claims."}
    (out / "stage2_checks.json").write_text(json.dumps(evidence, indent=2))
    print(json.dumps({"T6": white, "filtered": filtered, "T7_equal_arm_max": equal_arm_max}, indent=2))
    return 0 if all(r["pass"] for r in white) and equal_arm_max == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
