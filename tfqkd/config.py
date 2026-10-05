"""Single configuration entry point. Run: python -m tfqkd.config."""

import argparse
import json
from pathlib import Path
import tomllib
import warnings

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = ROOT / "configs/config.toml"


def load(path=DEFAULT_CONFIG):
    with Path(path).open("rb") as stream:
        config = tomllib.load(stream)
    for section in config.values():
        if isinstance(section, dict):
            for key, value in section.items():
                if key.startswith("PLACEHOLDER_") or str(value).startswith("PLACEHOLDER_"):
                    warnings.warn(f"Unresolved input: {key}={value}", RuntimeWarning)
    p, g, o = config["physics"], config["grid"], config["operation"]
    if not 0 < g["f_min_hz"] < g["f_max_hz"] or g["points"] < 2:
        raise ValueError("Frequency bounds must be positive and increasing; at least two points required")
    if p["n"] <= 0 or p["K"] not in (2, 4):
        raise ValueError("n must be positive; K must be 2 or 4 (Eq. 5 discussion)")
    if o["sigma_limit_rad"] <= 0 or o["tau_max_s"] <= 0 or o["tau_ps_s"] < 0:
        raise ValueError("Invalid operating threshold or times")
    if g.get('mode', 'reference') not in ('fast', 'reference'):
        raise ValueError('Grid mode must be fast or reference')
    if any(g.get(key, g['points']) < 2 for key in ('fast_points', 'reference_points')):
        raise ValueError('Each grid mode needs at least two points')
    if 'performance' in config:
        settings = config['performance']
        if settings['workers'] < 1 or settings['cache_entries'] < 1 or settings['batch_memory_mb'] <= 0:
            raise ValueError('Invalid process/cache/working-array configuration')
    if 'classical_actuator' in config:
        ac = config['classical_actuator']
        if not 0 < ac['omega_a_min_rad_s'] <= ac['omega_a_max_rad_s'] or ac['omega_a_points'] < 2:
            raise ValueError('Invalid engineering actuator-pole scan')
        if any(not 0 < fraction < 1 for fraction in ac['gain_fractions']):
            raise ValueError('Gain fractions must be strictly inside the stable interval (0,1)')
        if not 0 < ac['gain_ceiling_fraction'] < 1:
            raise ValueError('The gain scan must stop strictly below marginal stability')
        if not 0 < ac['length_min_km'] < ac['length_max_km'] or ac['length_points'] < 2:
            raise ValueError('Invalid balanced-arm length scan')
    return config


def cli_config():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument('--grid-mode', choices=('fast', 'reference'))
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    config = load(args.config)
    if args.grid_mode is not None:
        config['grid']['mode'] = args.grid_mode
    if args.output_dir is not None:
        config['_output_dir'] = str(args.output_dir.resolve())
    return config


if __name__ == "__main__":
    print(json.dumps(cli_config(), indent=2))
