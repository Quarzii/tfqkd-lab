"""Laser loop only; delayed fiber compensation belongs to stage 2."""

import numpy as np

from .config import cli_config


def positive_frequencies(f):
    values = np.asarray(f, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("Frequencies must be finite and strictly positive")
    return values


def laser_gain(f, p):
    f = positive_frequencies(f)
    # Eq. F4, bertaina2024; following paragraph; QKD.ipynb G1 (unrounded G0).
    G0 = (2 * np.pi * p["B_hz"]) ** 2 * (1 + p["delta"]) / (1 + p["gamma"])
    # Eq. F4, bertaina2024.
    return G0 / (2j * np.pi * f) ** 2 * (1j * f + p["B_hz"] * p["gamma"]) / (1j * f + p["B_hz"] * p["delta"])


def laser_residual_power(f, p):
    # Eq. F2, bertaina2024.
    return np.abs(1 / (1 + laser_gain(f, p))) ** 2


def common_laser_power(f, delta_L_km, p):
    f = positive_frequencies(f)
    # Eq. 5, bertaina2024; f [Hz], delta_L [km], c [km/s].
    return 4 * np.sin(2 * np.pi * f * p["n"] * delta_L_km / p["c_km_s"]) ** 2


if __name__ == "__main__":
    c = cli_config()
    f = np.geomspace(c["grid"]["f_min_hz"], c["grid"]["f_max_hz"], c["grid"]["points"])
    print("laser residual power min/max:", laser_residual_power(f, c["physics"]).min(), laser_residual_power(f, c["physics"]).max())
