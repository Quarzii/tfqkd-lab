# Part C: reports, sensitivity and installation comparisons

The validated physics and part A numerical core are unchanged. Existing four-rate
ceilings and conditional inverse requirements are included directly, with no new
key-rate or noise model.

## One installation

```bash
python -m tfqkd.lab_run examples/bertaina2024_table3.toml --output results/my_report
python -m tfqkd.lab_run examples/missing_r3.toml --output results/my_requirement
```

Outputs: `result.json`, `report.md`, standalone `report.html`, `spectrum.png`,
`cumulative_variance.png`. HTML embeds the scientific plots as SVG and needs
no service, script, font download or Markdown package. Compute and rendering
times are recorded separately. JSON includes the complete spectrum and
accumulated integrals used by the contribution report.

Reports follow this order:

1. Four rates, H and realized classical/dual ceiling fractions.
2. Selected key rate, tau_Q, duty and the existing input-range envelope.
3. Ranked signed sensitivity to individual upgrades.
4. Existing conditional equipment requirements for unmeasured amplitudes.
5. Optimized classical g, g_crit and gain-distance margin; manual selected g is separate.
6. Laser/fiber/detection contributions at the **same selected tau_Q**, spectrum
   and accumulated-variance plots.
7. Applicability, **Tool limitations only**, provenance and run notices.

Classical phase-detection noise of the round-trip beat measurement is absent.
`CLASSICAL_DETECTION_NOISE` identifies R_classical and its realized fraction as
optimistic and the comparison as biased in favor of classical compensation.
No new detector coefficient or noise floor is inserted.

`UNKNOWNS.md` now separates **Research-stage assumptions**, **Tool limitations**
and **Validation limits**. Equal-arm research geometry, independent stabilized
research lasers and the old actuator scans are not imposed on user inputs.
Publication-specific T4/B6 RMS discrepancies stay in validation reports; they
are not printed as the Tool limitations section of a user run.

## What “improved by two” means

These transformations were explicitly accepted by the user as engineering
what-if choices, not publication-derived equipment upgrades:

| Quantity | Applied change |
| --- | --- |
| Noise amplitudes r3, r2, C2–C4, l, s0 | Half the coefficient |
| Detector dark-count rate and error | Half the value |
| Line attenuation [dB/km] | Half the value; this is a hypothetical installation change, not a noise reduction |
| Signed arm imbalance | Half the magnitude, preserving the physical total length |
| Detector efficiency eta | Halve the shortfall: eta → 1 − (1 − eta)/2; e.g. 0.9 → 0.95 |
| Pulse frequency | Double; the gain is approximately proportional clock scaling, not noise suppression |
| Classical actuator omega_a | Double the pole frequency; the existing stability/g search is recalculated |

Cutoffs, wavelengths, n, protocol settings and operation thresholds/times have
no approved improvement direction and are not ranked. Inactive terms are not
ranked. A supplied explicit zero has no factor-two change. **Unmeasured amplitudes
are never ranked as zero-noise sensitivity**: their inverse requirement is shown
instead. No effect on inactive parameters is represented as a practical upgrade.

Each row changes one parameter only and records the actual before/after value,
units, new key rate, signed absolute gain and relative gain. Gains are sorted by
absolute gain; any worsening remains negative. A supplied Lorentz linewidth is
displayed in Hz while r2 is changed through the established conversion.

Automatic g is re-optimized for each change. If the main input explicitly fixes
manual g, sensitivity retains that g; the recommendation remains the separate
optimized classical comparison. The margin 1−g/g_crit is not a phase margin in
degrees or a guarantee for the real loaded actuator. Equal arms reuse identical
Williams ratios numerically, without changing their equation or the arm sum.

For range inputs, the existing endpoint/corner envelope and zero-imbalance
extremum are retained. The sensitivity ranking and headline H refer to the
explicit point values; they do not claim to optimize or bound those quantities
through an unexamined nonmonotone continuum.

## Two or three installations

Each installation is a complete explicit configuration. No preset or shared
hidden apparatus values are selected:

```toml
[comparison]
configurations = ["cavity_classical.toml", "cavity_dual.toml"]
labels = ["Classical round-trip", "Dual-band"]
working_total_lengths_km = [50.0, 100.0, 150.0, 200.0, 250.0]
workers = 2
```

Paths are relative to the comparison file. Exactly two or three configurations
with distinct labels are accepted. The working coordinates above are the
user-requested C2 coordinates, not physical equipment defaults.

```bash
python -m tfqkd.lab_run examples/part_c/compare_two.toml --output results/two
python -m tfqkd.lab_compare examples/part_c/compare_three.toml --output results/three
```

The comparison report includes configured-point rates/H, top sensitivity,
recommended classical gain and margin, dominant noise and missing-input status.
Each variant has its own complete report in `variant_1/`, etc.

The working-length table reports key rate, tau_Q, duty, and their ratios to the
first variant. Length means **physical total length**, with each variant's
signed imbalance kept fixed. Zero key denominators are infinite/undefined
as appropriate. Classical rows carry optimistic-estimate labels.

The working-length table preserves a supplied noise/imbalance range envelope
in JSON. A range over total/individual arm length has no unique meaning during
a fixed-total comparison; it is rejected with instructions to use separately
specified variants. It is not silently discarded.

Independent variants and working coordinates use processes. Within a variant,
gain calculations remain vectorized; nested process pools are avoided. Use
`workers=1` for environments that cannot create worker processes. Python API
callers creating pools need the standard `if __name__ == '__main__'` guard.

## Actuator CSV phase gate correction

`[fit].maximum_rms_phase_deg` defaults to **5.0 degrees**, explicitly approved
as an engineering fit criterion, not a value from a source. Override it globally
or with `actuator.response.maximum_rms_phase_deg`. The magnitude/PSD gate remains
0.5 dex. Output records the effective phase gate and whether the response metadata
overrode it. The two gates do not establish apparatus calibration uncertainty.

## Reproduction and validation

```bash
python scripts/validate/validate_stage3_c.py --quick
python scripts/validate/validate_stage3_c.py --outputs --regression
python scripts/validate/validate_stage3_c.py --full
```

Results are in `results/stage3_c/`; the full suite keeps the accepted T4 discrepancy
as a validation limit, never as spectral agreement. The reference index and
source links remain documented in [UNIVERSAL_INPUTS.md](UNIVERSAL_INPUTS.md).
This part implements scientific reports and comparisons; the local browser
interface and final English user guide remain part D.
