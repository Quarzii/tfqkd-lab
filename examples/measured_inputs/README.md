# Measured-input examples and conditional experiment comparisons

The eight Pittaluga/Zhou TOMLs execute explicit phase **brackets**, not
same-window phase measurements at each key-rate length. Each file records its
source and assumptions. Complete arm losses and both detector specifications
are explicit. Missing intrinsic error and, for Pittaluga, missing fEC use labelled
ideal limits. No Bertaina apparatus coefficient is substituted.

Pittaluga's 0.104 rad is from Fig.2f at605 km; CAL is at368.702 km.
Zhou's maximum0.099 rad is from Fig.3e–g on lambda_c at615.6 km, not from
the three different key-rate lengths. The assumed phase/length inequality,
reference/quantum equivalence for Zhou, and frame/duty idealization are
comparison-specific. These are not built-in apparatus defaults.

`regression/*_laser.csv` and `regression/*_line.csv` are **synthetic source-model
fixtures**, generated with the Table III equations on65537 points.
They are not digitized experiments. `direct_table3.toml` uses one such F1 laser
fixture and explicit Table III arm/protocol values to demonstrate direct mode.

```bash
python -m tfqkd.lab_run examples/measured_inputs/direct_table3.toml --output results/my_direct_run
python -m tfqkd.lab_run examples/measured_inputs/zhou2023_518.16_phase_bound.toml --output results/my_phase_run
```

See [validation](../../docs/VALIDATION.md) for conditional comparisons and
[measured inputs](../../docs/MEASURED_INPUTS.md) for configuration details.
