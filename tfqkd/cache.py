"""Bounded process-local cache of unchanged spectra on an exact frequency grid."""
from collections import OrderedDict
from functools import wraps
from hashlib import blake2b
import numpy as np
from .config import load

_CACHE = OrderedDict()
_HITS = 0
_MISSES = 0
_LIMIT = load()['performance']['cache_entries']


def configure_cache(entries):
    global _LIMIT
    if entries < 1:
        raise ValueError('Cache capacity must be positive')
    _LIMIT = int(entries)
    while len(_CACHE) > _LIMIT:
        _CACHE.popitem(last=False)


def clear_spectral_cache():
    global _HITS, _MISSES
    _CACHE.clear(); _HITS = 0; _MISSES = 0


def cache_info():
    return dict(entries=len(_CACHE), hits=_HITS, misses=_MISSES)


def cached_spectrum(function):
    @wraps(function)
    def wrapped(frequency, physics):
        global _HITS, _MISSES
        f = np.ascontiguousarray(frequency, dtype=float)
        key = (function.__name__, f.shape, blake2b(f.tobytes()).digest(), tuple(sorted(physics.items())))
        if key in _CACHE:
            _HITS += 1
            _CACHE.move_to_end(key)
            return _CACHE[key]
        _MISSES += 1
        result = np.asarray(function(f, physics), dtype=float)
        result.setflags(write=False)
        _CACHE[key] = result
        while len(_CACHE) > _LIMIT:
            _CACHE.popitem(last=False)
        return result
    return wrapped


if __name__ == '__main__':
    from .config import cli_config
    config = cli_config()
    configure_cache(config['performance']['cache_entries'])
    print('Process-local spectrum cache:', cache_info(), 'capacity:', _LIMIT)
