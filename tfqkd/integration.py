"""Finite-band Eq. 4 integration with an explicit lower cutoff."""

import numpy as np
from scipy.integrate import cumulative_trapezoid

from .config import cli_config


def grid_points(config, mode=None):
    g = config['grid']
    selected = g.get('mode', 'reference') if mode is None else mode
    if selected not in ('fast', 'reference'):
        raise ValueError('Grid mode must be fast or reference')
    return g.get('fast_points', g['points']) if selected == 'fast' else g.get('reference_points', g['points'])


def frequency_grid(config, points=None, mode=None):
    g = config["grid"]
    return np.geomspace(g["f_min_hz"], g["f_max_hz"], grid_points(config, mode) if points is None else points)


class PhaseIntegral:
    def __init__(self, frequency_hz, psd):
        self.f = np.asarray(frequency_hz, dtype=float)
        self.psd = np.asarray(psd, dtype=float)
        if self.psd.ndim < 1 or self.f.shape[-1] != self.psd.shape[-1] or self.f.shape[-1] < 2:
            raise ValueError("The last frequency and PSD axes must match")
        if self.f.ndim > 1:
            self.f = np.broadcast_to(self.f, self.psd.shape)
        if np.any(~np.isfinite(self.f)) or np.any(self.f <= 0) or np.any(np.diff(self.f, axis=-1) < 0) or np.any(self.f[..., -1] <= self.f[..., 0]):
            raise ValueError("Frequency grid must be finite, positive and increasing")
        if np.any(~np.isfinite(self.psd)) or np.any(self.psd < 0):
            raise ValueError("PSD must be finite and nonnegative")
        # Eq. 4, bertaina2024: integrate each spectrum ONCE, from high to low frequency.
        self.cumulative = -cumulative_trapezoid(self.psd[..., ::-1], self.f[..., ::-1], axis=-1, initial=0)[..., ::-1]

    def _locate(self, values, target, descending=False):
        """Vectorized binary interpolation lookup, including padded duplicate grid nodes."""
        if values.ndim == 1:
            return np.searchsorted(-values if descending else values,
                                   -target if descending else target, side='right') - 1
        left = np.zeros(self.psd.shape[:-1], dtype=int)
        right = np.full(self.psd.shape[:-1], self.f.shape[-1] - 1, dtype=int)
        while np.any(right - left > 1):
            middle = (left + right) // 2
            inside = self._take(values, middle) >= target if descending else self._take(values, middle) <= target
            left = np.where(inside, middle, left)
            right = np.where(inside, right, middle)
        return left

    def _take(self, values, index):
        if values.ndim == 1:
            return values[index]
        return np.take_along_axis(values, np.asarray(index)[..., None], axis=-1)[..., 0]

    def above(self, lower_hz):
        lower = np.asarray(lower_hz, dtype=float)
        if np.any(~np.isfinite(lower)) or np.any(lower < self.f[..., 0]):
            raise ValueError("Lower cutoff lies below frequency grid; extend the grid")
        clipped = np.minimum(lower, self.f[..., -1])
        if self.psd.ndim > 1:
            clipped = np.broadcast_to(clipped, self.psd.shape[:-1])
        index = np.clip(self._locate(self.f, clipped), 0, self.f.shape[-1] - 2)
        # Eq. 4, bertaina2024, linear interpolation inside the boundary trapezoid.
        left, right = self._take(self.psd, index), self._take(self.psd, index + 1)
        lo, hi = self._take(self.f, index), self._take(self.f, index + 1)
        fraction = np.divide(clipped - lo, hi - lo, out=np.zeros_like(clipped, dtype=float), where=hi > lo)
        at_lower = left + (right - left) * fraction
        # Eq. 4, bertaina2024, partial-bin trapezoid plus all complete higher bins.
        return self._take(self.cumulative, index + 1) + (at_lower + right) * (hi - clipped) / 2

    def variance(self, tau_s):
        tau = np.asarray(tau_s, dtype=float)
        if np.any(~np.isfinite(tau)) or np.any(tau <= 0):
            raise ValueError("Integration time must be finite and positive")
        # Eq. 4, bertaina2024.
        return self.above(1 / tau)

    def operating_time(self, sigma_limit_rad, tau_max_s):
        if sigma_limit_rad <= 0 or tau_max_s <= 0:
            raise ValueError("Threshold and maximum time must be positive")
        # Eq. 4 and Sec. V, bertaina2024: threshold expressed as variance.
        limit = sigma_limit_rad**2
        capped = self.variance(tau_max_s) <= limit
        if np.all(capped):
            if self.psd.ndim == 1:
                return float(tau_max_s), 'capped'
            return np.full(self.psd.shape[:-1], tau_max_s), np.full(self.psd.shape[:-1], 'capped')
        # Eq. 4, bertaina2024: locate the threshold in the monotone accumulated integral.
        # Binary lookup below operates on indices only; there is no iterative integration or root solver.
        index = np.clip(self._locate(self.cumulative, limit, descending=True), 0, self.f.shape[-1] - 2)
        # Eq. 4, bertaina2024: invert the partial-bin integral of the linearly interpolated PSD.
        # This is exact inverse interpolation of the same trapezoids, not a new integration rule.
        rem = np.maximum(limit - self._take(self.cumulative, index + 1), 0)
        p_left, p_right = self._take(self.psd, index), self._take(self.psd, index + 1)
        lo, hi = self._take(self.f, index), self._take(self.f, index + 1)
        slope = np.divide(p_right - p_left, hi - lo, out=np.zeros_like(p_left), where=hi > lo)
        # Eq. 4 numerical primitive: rem=p_right*h-slope*h^2/2, solved in a cancellation-safe form.
        denominator = p_right + np.sqrt(np.maximum(p_right**2 - 2 * slope * rem, 0))
        distance = np.divide(2 * rem, denominator, out=np.zeros_like(rem, dtype=float), where=denominator > 0)
        cutoff = hi - distance
        # Eq. 4 / Sec. V, bertaina2024: tau=1/f_lower, with the configured upper time cap.
        tau = np.where(capped, tau_max_s, np.minimum(1 / cutoff, tau_max_s))
        status = np.where(capped, 'capped', 'threshold')
        return (float(tau), str(status)) if self.psd.ndim == 1 else (tau, status)


if __name__ == "__main__":
    from .spectra import components
    c = cli_config()
    f = frequency_grid(c)
    for s in c["scenarios"]:
        integral = PhaseIntegral(f, components(f, s, c)["total"])
        tau, status = integral.operating_time(c["operation"]["sigma_limit_rad"], c["operation"]["tau_max_s"])
        print(s["name"], tau, "s", float(integral.variance(tau)), "rad^2", status)
