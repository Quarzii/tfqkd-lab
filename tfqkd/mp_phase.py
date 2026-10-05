"""Mode-pairing phase increments only; no key-rate formula or tracking model.

One-sided phase PSDs [rad^2/Hz], frequencies [Hz], intervals [s], lengths [km].
The unchanged TF spectrum functions supply F1--F4/Eq. 6 where requested.
"""
from dataclasses import dataclass
from pathlib import Path
import argparse
import json
import tomllib

import numpy as np
from scipy.integrate import trapezoid

from .spectra import free_fiber
from .transfers import laser_residual_power


@dataclass(frozen=True)
class PhaseComponent:
    name: str
    psd: object
    low_frequency_power: float
    f2_psd_limit: float | None
    source: str

    @property
    def needs_low_cutoff(self):
        # Eq. 1, didomenico2010, difference kernel: integrand ~ f^(2-power).
        return self.low_frequency_power >= 3


@dataclass(frozen=True)
class IntegrationGrid:
    positive_floor_hz: float
    upper_hz: float
    points: int
    interval_batch: int

    def frequencies(self):
        if not 0 < self.positive_floor_hz < self.upper_hz:
            raise ValueError("The numerical positive floor and upper frequency must be positive and ordered")
        if self.points < 2 or self.interval_batch < 1:
            raise ValueError("At least two frequency nodes and a positive interval batch are required")
        # Numerical quadrature nodes for the increment integral, Eq. 1, didomenico2010.
        return np.geomspace(self.positive_floor_hz, self.upper_hz, self.points)


def _positive(name, value, *, zero=False):
    value = float(value)
    if not np.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")
    return value


def laser_components(name, settings):
    """Explicit models only: free F1, stabilized F2--F4, or Lorentzian white frequency."""
    model = settings["model"]
    if model == "lorentzian_white_frequency":
        linewidth = _positive(name + ".linewidth_hz", settings["linewidth_hz"], zero=True)
        # Eq. 1 and white-frequency paragraph after Eq. 2, didomenico2010: FWHM=pi*h0.
        coefficient = linewidth / np.pi
        if not coefficient:
            return []
        return [PhaseComponent(name + ".white_frequency", lambda f: coefficient / f**2,
                               2, coefficient, "Eq. 1, didomenico2010; Eq. 3, zhang2025")]
    if model not in ("free", "stabilized"):
        raise ValueError(f"Unknown laser model {model!r}; no apparatus preset is supplied")
    r3 = _positive(name + ".r3", settings["r3"], zero=True)
    r2 = _positive(name + ".r2", settings["r2"], zero=True)
    fc = _positive(name + ".fc_hz", settings["fc_hz"])
    terms = []
    if model == "free":
        if r3:
            # Eq. F1, bertaina2024: isolate the divergent term rather than cutting all F1.
            terms.append(PhaseComponent(name + ".r3", lambda f: r3 / f**3,
                                        3, None, "Eq. F1, bertaina2024"))
        if r2:
            # Eq. F1, bertaina2024: this term has a convergent phase-increment integral.
            terms.append(PhaseComponent(name + ".r2", lambda f: r2 / f**2 * (fc / (f + fc))**2,
                                        2, r2, "Eq. F1, bertaina2024"))
        return terms
    # Eq. F2/F4, bertaina2024: use the existing controller without editing its implementation.
    controller = {key: settings[key] for key in ("B_hz", "gamma", "delta")}
    for key, value in controller.items():
        _positive(name + "." + key, value)
    if r3 or r2:
        # Eq. F1/F2, bertaina2024; |1/(1+G)|^2 is O(f^4) at zero by Eq. F4.
        def residual(f):
            return laser_residual_power(f, controller) * (r3 / f**3 + r2 / f**2 * (fc / (f + fc))**2)
        terms.append(PhaseComponent(name + ".servo_residual", residual, -1 if r3 else -2,
                                    0, "Eq. F1/F2/F4, bertaina2024"))
    for key, power in (("C4", 4), ("C3", 3), ("C2", 2)):
        amplitude = _positive(name + "." + key, settings[key], zero=True)
        if amplitude:
            # Eq. F3, bertaina2024: each cavity term keeps its own lower limit.
            terms.append(PhaseComponent(name + "." + key,
                                        lambda f, a=amplitude, exponent=power: a / f**exponent,
                                        power, amplitude if power == 2 else None,
                                        "Eq. F3, bertaina2024"))
    return terms


def fiber_component(name, settings):
    length = _positive(name + ".length_km", settings["length_km"], zero=True)
    amplitude = _positive(name + ".l", settings["l"], zero=True)
    cutoff = _positive(name + ".fc1_hz", settings["fc1_hz"])
    # Eq. 6, bertaina2024: a single transit, no Eq. 5 round-trip multiplier K.
    coefficient = amplitude * length
    return PhaseComponent(name, lambda f: free_fiber(f, length, {"l": amplitude, "fc1_hz": cutoff}),
                          2, coefficient, "Eq. 6, bertaina2024; Fig. 2, zeng2022")


def components_from_config(config):
    """Independent lasers and independent one-pass arm noise, as explicitly requested."""
    # Eq. 7, bertaina2024; Fig. 2/Box 1, zeng2022: two users transmit to Charlie once.
    return (laser_components("laser_A", config["laser_A"])
            + laser_components("laser_B", config["laser_B"])
            + [fiber_component("fiber_A", config["fiber_A"]),
               fiber_component("fiber_B", config["fiber_B"])])


def component_variance(component, delta_t_s, grid, t_track_s=None):
    intervals = np.asarray(delta_t_s, dtype=float)
    if np.any(~np.isfinite(intervals)) or np.any(intervals < 0):
        raise ValueError("Pair intervals must be finite and nonnegative")
    frequencies = grid.frequencies()
    if t_track_s is not None:
        t_track_s = _positive("T_track [s]", t_track_s)
    if component.needs_low_cutoff:
        if t_track_s is None:
            raise ValueError(f"{component.name}: phase-increment integral diverges at f -> 0; "
                             "explicit T_track is required to define f_min=1/T_track. "
                             "This is a band restriction, not a frequency-tracking transfer function.")
        # User-approved band convention, 2026-10-03; Eq. 1, didomenico2010 is integrated over this band.
        lower = 1 / t_track_s
        if lower >= grid.upper_hz:
            raise ValueError("1/T_track is at or above the upper frequency of the integration band")
        # Numerical quadrature of Eq. 1, didomenico2010: cover the entire explicitly requested band.
        # A cutoff below the grid's usual positive floor must not leave one large trapezoid.
        if lower < grid.positive_floor_hz:
            frequencies = np.geomspace(lower, grid.upper_hz, grid.points)
        else:
            frequencies = np.r_[lower, frequencies[frequencies > lower]]
    else:
        lower = 0.0
        if component.f2_psd_limit is None or not np.isfinite(component.f2_psd_limit):
            raise ValueError(f"{component.name}: a finite limit of f^2*S(f) at zero must be supplied")
    psd = np.asarray(component.psd(frequencies), dtype=float)
    if psd.shape != frequencies.shape or np.any(~np.isfinite(psd)) or np.any(psd < 0):
        raise ValueError(f"{component.name}: PSD must be finite, nonnegative and match the frequency nodes")
    # Eq. 1, didomenico2010: rewrite 4*S*sin^2(pi*f*dt) to use the removable zero limit.
    f2_psd = frequencies**2 * psd
    nodes = frequencies if lower else np.r_[0.0, frequencies]
    if not lower:
        # Eq. 1, didomenico2010; Eq. 6/F1/F3, bertaina2024: exact f^2*S limit of each model term.
        f2_psd = np.r_[component.f2_psd_limit, f2_psd]
    flat = intervals.ravel()
    values = np.empty_like(flat)
    for start in range(0, len(flat), grid.interval_batch):
        batch = flat[start:start + grid.interval_batch, None]
        # Eq. 1, didomenico2010, algebraic phase-difference filter; np.sinc(x)=sin(pi*x)/(pi*x).
        integrand = 4 * np.pi**2 * batch**2 * f2_psd * np.sinc(batch * nodes)**2
        # Eq. 1, didomenico2010: numerical integral, with a zero-frequency endpoint for convergent terms.
        values[start:start + len(batch)] = trapezoid(integrand, nodes, axis=-1)
    return values.reshape(intervals.shape), lower


def phase_difference_variance(components, delta_t_s, grid, t_track_s=None):
    components = tuple(components)
    if len({c.name for c in components}) != len(components):
        raise ValueError("Phase component names must be unique")
    values, lower = {}, {}
    for component in components:
        values[component.name], lower[component.name] = component_variance(component, delta_t_s, grid, t_track_s)
    # Eq. 7, bertaina2024; Eq. 1, didomenico2010: independent one-pass PSDs and their integrals add.
    total = sum(values.values(), np.zeros_like(np.asarray(delta_t_s, dtype=float)))
    return {"variance_rad2": total, "components_rad2": values, "lower_bounds_hz": lower,
            "cut_components": [c.name for c in components if c.needs_low_cutoff],
            "upper_hz": grid.upper_hz, "T_track_s": t_track_s}


def tracking_band_sensitivity(components, delta_t_s, grid, t_track_s):
    t_track_s = _positive("T_track [s]", t_track_s)
    # User-requested band-sensitivity factors, 2026-10-03; no estimator model is added.
    return {label: phase_difference_variance(components, delta_t_s, grid, t_track_s * factor)
            for label, factor in (("shorter_x10", 0.1), ("nominal", 1.0), ("longer_x10", 10.0))}


def fiber_short_interval(delta_t_s, l, length_km, fc1_hz, *, next_order=False):
    t = np.asarray(delta_t_s, dtype=float)
    if np.any(~np.isfinite(t)) or np.any(t < 0):
        raise ValueError("Pair intervals must be finite and nonnegative")
    for name, value in (("l", l), ("length_km", length_km), ("fc1_hz", fc1_hz)):
        _positive(name, value, zero=name != "fc1_hz")
    # Eq. 6, bertaina2024 + Eq. 1, didomenico2010: derived small-interval second moment.
    result = 4 * np.pi**2 * l * length_km * fc1_hz * t**2
    if next_order:
        # Same Eqs.; derived tail correction: D/D_leading=1-(pi^2/3)*fc1*|dt|+O((fc1*dt)^2*|log(fc1*dt)|).
        result = result * (1 - np.pi**2 / 3 * fc1_hz * t)
    return result


def equivalent_fiber_drift(l, length_km, fc1_hz):
    for name, value in (("l", l), ("length_km", length_km), ("fc1_hz", fc1_hz)):
        _positive(name, value, zero=name != "fc1_hz")
    # Eq. 6, bertaina2024 + Eq. 3, zhang2025: derived equivalent small-interval Gaussian drift rate.
    return 2 * np.pi * np.sqrt(l * length_km * fc1_hz)


def zhang_x_error(variance_rad2, phase_slices):
    variance = np.asarray(variance_rad2, dtype=float)
    if np.any(~np.isfinite(variance)) or np.any(variance < 0):
        raise ValueError("Phase variance must be finite and nonnegative")
    if not isinstance(phase_slices, (int, np.integer)) or phase_slices < 2:
        raise ValueError("At least two phase slices are required")
    # Eq. B10, zhang2025: raw coherent-state X error; not the single-photon privacy error.
    return 0.5 - phase_slices / (8 * np.pi) * np.sin(2 * np.pi / phase_slices) * np.exp(-variance / 2)


def zhou_x_error(variance_rad2, delta_t_s, visibility_v2, residual_frequency_hz):
    variance, interval = np.broadcast_arrays(np.asarray(variance_rad2, dtype=float),
                                            np.asarray(delta_t_s, dtype=float))
    if np.any(~np.isfinite(variance)) or np.any(variance < 0) or np.any(~np.isfinite(interval)) or np.any(interval < 0):
        raise ValueError("Phase variance and pair intervals must be finite and nonnegative")
    if not np.isfinite(visibility_v2) or not 0 <= visibility_v2 <= 0.5:
        raise ValueError("V2 must lie in [0,0.5] for the coherent-state model in Zhou Eq. 2")
    if not np.isfinite(residual_frequency_hz):
        raise ValueError("Residual beat frequency [Hz] must be specified and finite")
    # Eq. 2, zhou2023_async: deterministic residual beat is separate from the stochastic variance.
    return 0.5 - visibility_v2 / 2 * np.exp(-variance / 2) * np.cos(2 * np.pi * residual_frequency_hz * interval)


def discrete_average(error, weights):
    error, weights = np.asarray(error, dtype=float), np.asarray(weights, dtype=float)
    if error.shape != weights.shape or error.ndim != 1 or not len(error):
        raise ValueError("Errors and explicit interval counts/weights must be matching nonempty vectors")
    if np.any(~np.isfinite(error)) or np.any(~np.isfinite(weights)) or np.any(weights < 0) or not np.any(weights > 0):
        raise ValueError("Errors/weights must be finite; nonnegative weights must have a positive sum")
    if np.any(error < 0) or np.any(error > 1):
        raise ValueError("X error probabilities must lie in [0,1]")
    # Eq. D5, zhang2025, following paragraph recommends the discrete sum; weights here are explicit counts.
    return np.sum(error * weights) / np.sum(weights)


def _serializable(value):
    if isinstance(value, dict):
        return {k: _serializable(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with args.config.open("rb") as stream:
        config = tomllib.load(stream)
    components = components_from_config(config)
    grid = IntegrationGrid(**config["grid"])
    intervals = config["operation"]["delta_t_s"]
    track = config["operation"].get("T_track_s")
    result = phase_difference_variance(components, intervals, grid, track)
    if track is not None:
        result["tracking_band_sensitivity"] = tracking_band_sensitivity(components, intervals, grid, track)
    result["scope"] = "Phase increments only; no key-rate formula; independent one-pass components"
    result["low_cutoff_convention"] = "1/T_track only for divergent components; no frequency-tracking model"
    rendered = json.dumps(_serializable(result), indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
