"""Kasdin FIR generation and its exact discrete spectral density."""

import numpy as np
from scipy.signal import fftconvolve
from scipy.special import gamma
from .config import cli_config


def coefficients(samples, alpha):
    if samples < 2 or not 0 < alpha < 1:
        raise ValueError("This validation uses the stationary case 0 < alpha < 1")
    h = np.ones(samples)
    k = np.arange(1, samples)
    # Eq. 104, kasdin1995: h0=1; hk=hk-1*(alpha/2+k-1)/k.
    h[1:] = np.cumprod((alpha / 2 + k - 1) / k)
    return h


def generate(rng, samples, alpha, input_variance):
    if input_variance <= 0:
        raise ValueError("Input variance must be positive")
    h = coefficients(samples, alpha)
    # Appendix II and Eq. 37, kasdin1995: Gaussian IID input, variance Qd.
    white = rng.normal(scale=np.sqrt(input_variance), size=samples)
    # Eq. 37, Sec. III-B and VI-C, Appendix II, kasdin1995: LINEAR convolution.
    return fftconvolve(h, white, mode="full")[:samples]


def discrete_psd_one_sided(f, fs, alpha, input_variance):
    f = np.asarray(f, dtype=float)
    if np.any(f <= 0) or np.any(f > fs / 2):
        raise ValueError("Require 0 < f <= Nyquist")
    # Eq. 98 with bilateral convention Eq. 15, kasdin1995; fold negative frequencies.
    psd = 2 * input_variance / fs / (2 * np.sin(np.pi * f / fs)) ** alpha
    # Eq. 98, kasdin1995, unique Nyquist bin is not doubled in discrete one-sided PSD.
    return np.where(f == fs / 2, psd / 2, psd)


def power_law_one_sided(f, fs, alpha, input_variance):
    # Eq. 99 and bilateral convention Eq. 15, kasdin1995; positive-frequency folding.
    return 2 * input_variance * (1 / fs) ** (1 - alpha) / (2 * np.pi * f) ** alpha


def stationary_variance(alpha, input_variance):
    # Eq. 111, kasdin1995, valid for alpha < 1.
    return input_variance * gamma(1 - alpha) / gamma(1 - alpha / 2) ** 2


def finite_record_expectations(samples, alpha, input_variance):
    h = coefficients(samples, alpha)
    # Eq. 38 at lag m=0, kasdin1995; average the finite-start variances over time.
    mean_square = input_variance * np.mean(np.cumsum(h**2))
    # Eqs. 37-38, kasdin1995; derived variance of the record mean from IID input weights.
    mean_variance = input_variance * np.sum(np.cumsum(h)**2) / samples**2
    # Eqs. 37-38, kasdin1995; E[mean((x-mean(x))^2)] = E[mean(x^2)] - E[mean(x)^2].
    return {"mean_square": float(mean_square), "centered_variance": float(mean_square - mean_variance),
            "variance_of_record_mean": float(mean_variance)}


if __name__ == "__main__":
    v = cli_config()["validation"]
    x = generate(np.random.default_rng(v["seed"]), v["samples"], v["alpha"], v["input_variance"])
    print("sample variance:", np.var(x), "stationary Eq.111:", stationary_variance(v["alpha"], v["input_variance"]))
