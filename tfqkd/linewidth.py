"""Single-sided frequency PSD to optical linewidth: local Di Domenico source."""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import sici

from .config import cli_config


def white_linewidth(h0):
    if not np.isfinite(h0) or h0 <= 0:
        raise ValueError("Require a positive single-sided white frequency PSD in Hz²/Hz")
    # Eq. 5 and following text, didomenico2010, p. 4803: FWHM = pi*h0.
    return np.pi * h0


def normalized_autocorrelation(u, cutoff_ratio=None):
    # Eqs. 1, 4, 5, didomenico2010: u=h0*|tau|, optical carrier removed.
    if cutoff_ratio is None:
        return np.exp(-np.pi**2 * np.abs(u))
    # Eq. 4, didomenico2010; fc/h0 is dimensionless; E0 cancels in normalized spectra.
    argument = np.pi * cutoff_ratio * u
    return np.exp(2 / cutoff_ratio * (np.sin(argument)**2 - argument * sici(2 * argument)[0]))


def numerical_linewidth(h0, tolerance, cutoff_hz=None):
    white_linewidth(h0)
    if cutoff_hz is not None and (not np.isfinite(cutoff_hz) or cutoff_hz <= 0):
        raise ValueError("Cutoff must be positive")
    # Eqs. 2, 4, didomenico2010; change variables u=h0*tau and x=detuning/h0.
    ratio = None if cutoff_hz is None else cutoff_hz / h0

    def line(x):
        # Eq. 2, didomenico2010: even autocorrelation, positive-lag cosine integral.
        value, error = quad(lambda u: normalized_autocorrelation(u, ratio), 0, np.inf,
                            weight="cos", wvar=2 * np.pi * x, epsabs=tolerance)
        return value

    center = line(0)
    # Eqs. 2 and 5, didomenico2010: numerical half-maximum root; bracket from white limit.
    half_width = brentq(lambda x: line(x) - center / 2, 0, np.pi, xtol=tolerance)
    return 2 * half_width * h0


if __name__ == "__main__":
    config = cli_config()["linewidth"]
    for h0 in config["white_levels_hz2_per_hz"]:
        print("h0 [Hz²/Hz]", h0, "FWHM analytic/numerical [Hz]",
              white_linewidth(h0), numerical_linewidth(h0, config["quadrature_tolerance"]))
