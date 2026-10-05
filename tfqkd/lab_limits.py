"""User-facing limitations come ONLY from Tool limitations in UNKNOWNS.md."""
from .config import ROOT

CLASSICAL_DETECTION_NOISE = ('CLASSICAL_DETECTION_NOISE: phase-detection noise of the round-trip beat measurement is not included; '
    'the supplied sources do not parameterize it. R_classical and its realized ceiling fraction are optimistic estimates; '
    'comparison with R_dual is biased in favor of the classical scheme.')

def tool_limitations_markdown():
    text=(ROOT/'docs/UNKNOWNS.md').read_text()
    return text.split('## Tool limitations\n',1)[1].split('\n## ',1)[0].strip()

def limitation_ids():
    return [line.split('|')[1].strip() for line in tool_limitations_markdown().splitlines()
            if line.startswith('| ') and not line.startswith('| Identifier')]
