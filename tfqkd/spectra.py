"""One-sided phase PSDs, rad^2/Hz; lengths in km."""

import numpy as np

from .config import cli_config
from .transfers import common_laser_power, laser_residual_power, positive_frequencies
from .cache import cached_spectrum


@cached_spectrum
def free_laser(f, p):
    f = positive_frequencies(f)
    # Eq. F1, bertaina2024.
    return p["r3"] / f**3 + p["r2"] / f**2 * (p["fc_hz"] / (f + p["fc_hz"])) ** 2


@cached_spectrum
def cavity(f, p):
    f = positive_frequencies(f)
    # Eq. F3, bertaina2024.
    return p["C4"] / f**4 + p["C3"] / f**3 + p["C2"] / f**2


@cached_spectrum
def stabilized_laser(f, p):
    # Eq. F2, bertaina2024.
    return cavity(f, p) + laser_residual_power(f, p) * free_laser(f, p)


def free_fiber(f, length_km, p):
    f = positive_frequencies(f)
    length = np.asarray(length_km, dtype=float)
    if np.any(~np.isfinite(length)) or np.any(length < 0):
        raise ValueError("Fiber length must be finite and nonnegative")
    # Eq. 6, bertaina2024.
    factor = length if length.ndim == 0 else length[..., None]
    return factor * free_fiber_per_km(f, p)


@cached_spectrum
def free_fiber_per_km(f, p):
    # Eq. 6, bertaina2024: exact linear dependence on L; cached coefficient per km.
    return p['l'] / f**2 * (p['fc1_hz'] / (f + p['fc1_hz']))**2


def stabilized_fiber(f, length_km, p):
    f = positive_frequencies(f)
    length = np.asarray(length_km, dtype=float)
    if np.any(~np.isfinite(length)) or np.any(length < 0):
        raise ValueError("Fiber length must be finite and nonnegative")
    # Eq. 8, bertaina2024; QKD.ipynb stabfibnoise, excludes detection floor.
    factor = length if length.ndim == 0 else length[..., None]
    return factor * stabilized_fiber_per_km(f, p)


@cached_spectrum
def stabilized_fiber_per_km(f, p):
    # Eq. 8, bertaina2024: exact length coefficient, excluding detection noise.
    return ((p['lambda_s_nm'] - p['lambda_q_nm']) / p['lambda_s_nm'])**2 * p['l'] / f**2


@cached_spectrum
def detection(f, p):
    f = positive_frequencies(f)
    # Appendix G (unnumbered S_detection), bertaina2024; QKD.ipynb fiberdetection.
    return p["s0"] * (p["fc2_hz"] / (f + p["fc2_hz"])) ** 2


def components(f, scenario, config, LB_km=None):
    p = config["physics"]
    LB = config["operation"]["LB_km"] if LB_km is None else LB_km
    # Eq. 5, bertaina2024, definition of delta L; QKD.ipynb calc_spectra.
    LA = LB + scenario["delta_L_km"]
    laser_fn = stabilized_laser if scenario["cavity"] else free_laser
    fiber_fn = stabilized_fiber if scenario["stabilized"] else free_fiber
    if scenario["common"]:
        # Eq. 5, bertaina2024: K=4 correlated round trip, K=2 otherwise.
        laser = common_laser_power(f, scenario["delta_L_km"], p) * laser_fn(f, p)
        fiber = p["K"] * (fiber_fn(f, LA, p) + fiber_fn(f, LB, p))
    else:
        # Eq. 7, bertaina2024: no round trip, fiber coefficient is 1, not K.
        laser = laser_fn(f, p) + laser_fn(f, p)
        fiber = fiber_fn(f, LA, p) + fiber_fn(f, LB, p)
    detector = detection(f, p) if scenario["stabilized"] else np.zeros_like(f, dtype=float)
    # Eqs. 5, 7 and Appendix G, bertaina2024; QKD.ipynb calc_spectra: detector once.
    return {"laser": laser, "fiber": fiber, "detection": detector, "total": laser + fiber + detector}


if __name__ == "__main__":
    c = cli_config()
    from .integration import frequency_grid
    f = frequency_grid(c)
    for scenario in c["scenarios"]:
        s = components(f, scenario, c)["total"]
        print(scenario["name"], "PSD at grid endpoints [rad^2/Hz]:", s[0], s[-1])
