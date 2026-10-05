"""Current block 5: actuator pole, adaptive tau_Q and key-rate reach."""

# CLI import bootstrap; no calculation settings are changed.
import sys as _sys
from pathlib import Path as _Path
_PROJECT_ROOT = _Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_PROJECT_ROOT))
ROOT = _PROJECT_ROOT
from scripts.run.run_classical_actuator import main

if __name__ == "__main__":
    main()
