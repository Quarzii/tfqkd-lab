# Mode-pairing examples

Run from the repository root with dependencies installed.

## Phase increments

```bash
python -m tfqkd.mp_phase examples/mode_pairing/phase.toml --output results/my_mp_phase.json
```

`phase.toml` supplies explicit laser/fiber spectra, pairing intervals, and a
tracking-band convention. JSON reports phase-difference variance and tracking
window sensitivity. Bertaina fiber coefficients are transferred demonstration
inputs, not measured PSDs of the mode-pairing experiments. This command does
not calculate key rate.

## Asymptotic key rate

```bash
python -m tfqkd.mp_key examples/mode_pairing/key_zhang202_given_gains.toml --output results/my_mp_key
```

The example uses Zhang's observed gains at 202 km and explicit phase/protocol
parameters. It estimates a base-protocol key bound conditional on those gains;
it is not a forecast from fiber length and detector efficiency alone.

See the [user guide](../../docs/USER_GUIDE.md) for normalization and units,
[validation](../../docs/VALIDATION.md) for discrepancies, and
[source counts](../../sources/data/mode_pairing/key_comparison_inputs.json).
