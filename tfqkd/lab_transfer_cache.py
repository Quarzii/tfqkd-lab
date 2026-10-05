"""Exact reuse of Williams transfer ratios in laboratory sensitivity reports.

Only scheduling changes: keys include every frequency, delay and controller
parameter. No ratio is approximated. The existing batch-memory budget bounds
retained arrays; callers must treat cached results as immutable.
"""
from collections import OrderedDict
from hashlib import sha256
from threading import Lock
import numpy as np
from .classical import remote_controller_ratio

_VALUES = OrderedDict()
_LOCK = Lock()


def search_ratio(frequency, one_way_s, controller):
    """Same A8 arithmetic, with gain-independent spatial factors before broadcast.

    Used only for the search's one-dimensional log grid. The selected winner
    still goes through the original core quadrature and transfer function.
    """
    f = np.asarray(frequency, dtype=float)
    if f.ndim != 1 or np.any(f <= 0) or one_way_s < 0:
        raise ValueError('Search ratio requires a positive one-dimensional grid and nonnegative delay')
    # Eqs. A3,A6,A8, williams2008: identical expressions to remote_controller_ratio.
    s = 2j * np.pi * f
    C = controller(s)
    G = C * (1 + np.exp(-2 * s * one_way_s))
    u = 2 * np.pi * f * one_way_s
    D = 2 * C / (1 + G) * np.exp(-2 * s * one_way_s)
    # Eqs. A3,A4,A8, williams2008: spatial factors are independent of controller gain.
    cos2 = (1 + np.sinc(2 * u / np.pi)) / 2
    cross = (1 + np.exp(-1j * u) * np.sinc(u / np.pi)) / 2
    result = 1 + np.abs(D)**2 * cos2 - 2 * np.real(np.conj(D) * cross)
    small = np.abs(u) < 1
    if np.any(small):
        us = u[small]
        sin2 = np.zeros_like(us)
        # Eqs. A3,A8, williams2008: the core's stable Taylor evaluation, same stopping rule.
        term = us**2 / 3
        order = 1
        while np.any(np.abs(term) > np.finfo(float).eps * np.maximum(sin2, np.finfo(float).tiny)):
            sin2 += term
            # Eqs. A3,A8, williams2008: successive terms of the exact spatial integral.
            term *= -(2*us)**2 / ((2*order+2)*(2*order+3))
            order += 1
        # Eqs. A3,A8, williams2008: same low-frequency expression as the original core.
        b = 1-D[..., small]
        result[..., small] = np.abs(b)**2*(1-sin2)+sin2-np.imag(b)*us*np.sinc(us/np.pi)**2
    if np.any(result < 0) or np.any(~np.isfinite(result)):
        raise ArithmeticError('Invalid remote PSD ratio')
    return result


def clear():
    with _LOCK:
        _VALUES.clear()


def cached_ratio(frequency, one_way_s, controller, gains, omega_a, memory_mb, *, search=False):
    f = np.ascontiguousarray(frequency)
    g = np.ascontiguousarray(gains)
    key = (search, f.shape, f.dtype.str, sha256(f.tobytes()).digest(),
           float(one_way_s), g.shape, g.dtype.str, g.tobytes(), float(omega_a))
    budget = int(memory_mb * 1024**2)
    with _LOCK:
        if key in _VALUES:
            _VALUES.move_to_end(key)
            return _VALUES[key]
    # Eq.A8, williams2008: unchanged analytic spatially integrated remote ratio.
    value = (search_ratio(f, one_way_s, controller) if search else
             remote_controller_ratio(f, one_way_s, controller))
    value.setflags(write=False)
    if value.nbytes <= budget:
        with _LOCK:
            _VALUES[key] = value
            _VALUES.move_to_end(key)
            while sum(item.nbytes for item in _VALUES.values()) > budget:
                _VALUES.popitem(last=False)
    return value


if __name__ == '__main__':
    print('Exact transfer-ratio cache; frequency, delay, gain and pole are all part of the key.')
