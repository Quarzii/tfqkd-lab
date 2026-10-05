# Complete laboratory examples

All values are explicit calculation inputs, not equipment presets.

- `diode_classical.toml`: free diode laser, classical round-trip compensation.
- `stabilized_dual.toml`: cavity-stabilized independent lasers, dual-band compensation.
- `actuator_comparison_10km_arms.toml`: two 10 km arms, poles 2*pi*100 / 2*pi*1e5 rad/s; the fast pole gives a visible gain and g>0.
- `actuator_comparison_100km_arms.toml`: two 100 km arms, the same poles; both optimizers choose g=0 and the rates coincide.

The pole values and geometry are user-requested engineering examples, not measured specifications.
Laser/fiber coefficients: [Bertaina2024 Table III, Appendices F/G](https://doi.org/10.1002/qute.202400032).
Receiver/protocol: Table II, Appendix D and [QKD.ipynb](../../sources/data/QKD.ipynb), Cells 3/23.
Classical rates remain optimistic because phase-detection noise is omitted.

Run any configuration from the project root, for example:

```bash
python -m tfqkd.lab_run examples/laboratory/actuator_comparison_10km_arms.toml --output results/short_arms
python -m tfqkd.lab_run examples/laboratory/actuator_comparison_100km_arms.toml --output results/long_arms
python scripts/report/report_release_examples.py
```

The last command reproduces both committed reports with the checked example
interpretation under `reports/actuator_comparison_*km_arms/`.
At g=0 the Williams Eq.A8 remote residual is the uncompensated fiber PSD,
independent of the actuator pole. This explains equal outputs for the 100 km
arms: under these particular inputs the converged stable-gain search finds no
key improvement at its 1e-3 rate resolution. It is not a claim that every window
fails the phase threshold or that every installation has the same crossover.
Reach is off by default; the operating time is solved, not fixed at the 100 ms cap.
