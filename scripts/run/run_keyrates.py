"""Generate Stage 2 Figure 3 key-rate data."""

from __future__ import annotations

# CLI import bootstrap; no calculation settings are changed.
import sys as _sys
from pathlib import Path as _Path
_PROJECT_ROOT = _Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_PROJECT_ROOT))
ROOT = _PROJECT_ROOT

import argparse
import csv
import json
from pathlib import Path

import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/kvant-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tfqkd.config import load
from tfqkd.keyrates import (
    KeyRateParameters,
    config_toml_snippet,
    figure3_panel_scenarios,
    rates_for_loss,
    table_ii_sns_pd_parameters,
)


def positive_or_nan(values):
    array = np.asarray(values, dtype=float)
    return np.where(array > 0, array, np.nan)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(ROOT / "configs/config.toml"))
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--loss-min-db", type=float, default=15.0)
    parser.add_argument("--loss-max-db", type=float, default=110.0)
    parser.add_argument("--loss-points", type=int, default=200)
    args = parser.parse_args()

    config = load(args.config)
    params = KeyRateParameters(**config["keyrate"])
    args.loss_min_db = config["keyrate_validation"]["loss_min_db"]
    args.loss_max_db = config["keyrate_validation"]["loss_max_db"]
    args.loss_points = config["keyrate_validation"]["loss_points"]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    losses = np.linspace(args.loss_min_db, args.loss_max_db, args.loss_points)
    rows = []
    for panel, scenario in figure3_panel_scenarios(config).items():
        for loss_db in losses:
            row = rates_for_loss(config, scenario, float(loss_db), params=params)
            rows.append({"panel": panel, "scenario": scenario["name"], **row})

    csv_path = outdir / "keyrate_figure3.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axs = plt.subplots(2, 2, figsize=(10.5, 5.2), sharex=True, sharey=True)
    panel_axes = {"a": axs[0, 0], "b": axs[0, 1], "c": axs[1, 0], "d": axs[1, 1]}
    panel_titles = {"a": "a) Scenarios 1, 4", "b": "b) Scenarios 2, 5, 7", "c": "c) Scenario 3", "d": "d) Scenario 6"}
    for panel, ax in panel_axes.items():
        panel_rows = [row for row in rows if row["panel"] == panel]
        x = np.array([row["loss_db"] for row in panel_rows])
        # Figure 3, bertaina2024; QKD.ipynb Cell 23 plot_panel_scenarios.
        ax.semilogy(x, positive_or_nan([row["plob_bps"] for row in panel_rows]), label="Realistic PLOB bound", linewidth=1.5)
        # Figure 3, bertaina2024; QKD.ipynb Cell 23 plot_panel_scenarios.
        ax.semilogy(x, positive_or_nan([row["bb84_bps"] for row in panel_rows]), label="Decoy states BB84", linestyle=(0, (1, 1)), linewidth=2.0)
        # Figure 3, bertaina2024; QKD.ipynb Cell 23 plot_panel_scenarios.
        ax.semilogy(x, positive_or_nan([row["sns_aopp_bps"] for row in panel_rows]), label="SNS-AOPP", linestyle=(0, (4, 1)), linewidth=2.0)
        # Figure 3, bertaina2024; QKD.ipynb Cell 23 plot_panel_scenarios.
        ax.semilogy(x, positive_or_nan([row["cal_bps"] for row in panel_rows]), label="CAL", linestyle=(0, (1, 2, 4, 2)), linewidth=2.0)
        ax.set_xlim(args.loss_min_db, args.loss_max_db)
        ax.set_ylim(1e1, 1e7)
        ax.grid(color=(0.9, 0.9, 0.9))
        ax.tick_params(direction="in", which="both")
        ax.text(18, 50, panel_titles[panel], backgroundcolor="white")
    axs[1, 0].set_xlabel("Alice-Bob channel loss [dB]")
    axs[1, 1].set_xlabel("Alice-Bob channel loss [dB]")
    axs[0, 0].set_ylabel("Key rate [bit s$^{-1}$]")
    axs[1, 0].set_ylabel("Key rate [bit s$^{-1}$]")
    axs[0, 0].legend(fontsize=9, loc=(0.5, 0.5), frameon=True)
    for ax in axs[0]:
        ax.secondary_xaxis("top", functions=(lambda loss: loss / params.attenuation_db_per_km, lambda length: length * params.attenuation_db_per_km)).set_xlabel("Equivalent fiber length [km]")
    fig.tight_layout()
    png_path = outdir / "keyrate_figure3.png"
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    convergence = []
    for scenario in [figure3_panel_scenarios(config)["a"], figure3_panel_scenarios(config)["b"]]:
        for loss_db in (40.0, 80.0):
            for points in (4097, 16385, config["grid"]["points"]):
                row = rates_for_loss(config, scenario, loss_db, params=params, points=points)
                convergence.append({"scenario": scenario["name"], "loss_db": loss_db, "points": points, **row})

    json_path = outdir / "keyrate_evidence.json"
    json_path.write_text(json.dumps({"config_snippet": config_toml_snippet(params), "convergence": convergence}, indent=2))
    print(json.dumps({"csv": str(csv_path), "png": str(png_path), "json": str(json_path), "rows": len(rows)}, indent=2))


if __name__ == "__main__":
    main()
