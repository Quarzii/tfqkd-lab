# Mode-pairing key implementation and executed validation

Scope: base asymptotic Zeng 2022 protocol only. Separate new module `tfqkd/mp_key.py`; no TF kernel, accepted MP phase code, earlier tests or configuration changed (57 baseline SHA-256 matches). No async vacuum-key term or click-filter variant. No Zhu phase forecast.

## Implemented equations and required inputs

- Zeng 2022 Eq.4 / Supplement72: `r_p` from explicit single-click probability p and integer L_max; stable evaluation also for pL much smaller than one.
- Zeng 2022 Eq.7: r_p*r_s*[q11*(1-h(e11X))-f_EC*h(e_Z)], bits per shared quantum emission round. r_s, e_Z and f_EC are explicit.
- Zeng 2022 Eq.15/16: successful IID gaps use the conditional discrete geometric distribution. Zhang D5: count-weighted discrete error average; measured intervals may replace the IID law.
- Accepted MP increment integral: one-pass laser and fiber PSDs; no TF K. Spectral mode uses the unchanged phase module and explicit grid / T_track when a divergent term requires it. Published-drift mode uses Zhang Eq.3 (equal independent laser linewidths) and the explicit sigma_L.
- Zhang B10 and D2: raw coherent X error with finite phase slices and explicit intrinsic error. This is the weak balanced setting approximation, not a detector/dark-count model and not a single-photon privacy error.
- Zhang C4/C5/C6/C7: source asymptotic analytical decoy yield/error bounds; finite-key sampling theta_X=0. A nonpositive yield lower bound is vacuous; an error upper interval containing1/2 gives zero privacy. No25% subtraction. Negative error-yield upper bounds are rejected.
- Zeng Supplement54/68/69: q11 from the signal Poisson one-photon weight and its decoy yield lower bound. Gains/error gains are required matrices, not filled from a different apparatus.
- Zhu Eq.1 following paragraph: N_pair=N_rounds/2, so its published bit/potential-pair rate is twice the shared-round rate. Zhu TableVI and Zhang setup/AppendixE supply explicit clock/reference/recovery timing.

All physical expressions have source comments in code. Numerical fixtures/tolerances are test choices, not equipment coefficients.

## Input provenance and source conversion

Full original parameters, all counters, source locations and missing inputs: `sources/data/mode_pairing/key_comparison_inputs.json`. The calculation inputs and every exposure matrix are in `details.json`.

Zhang: TablesI–III, PDF p5,13–14, Eq.3, AppendixB/E. The source simulation uses linewidth 100 Hz for each laser and sigma_L=4000 rad/s; this is not a per-length measured drift spectrum. Active clock 1.25 GHz, reference 1/8. Single-click p is not tabulated: the source Zeng Supplement73/74 approximation is evaluated with the published receiver/loss data. Unfitted e_d=0 is explicitly ideal; calibration changes only D2 e_d.

Zhu: TablesI/VI/VII, PDF p4,16–17. Width2 kHz does not supply drift or tracking residuals. Per user mu/nu/probabilities and observed p/counts are tabulated. Quantum clock 444.83 MHz after strong 160.97 MHz and recovery 19.20 MHz out of625 MHz. No phase prediction is run, per user decision.

Source settings: Z per-user probabilities [P0²,2P0Pnu,2P0Pmu]; X [P0²,Pnu²,Pmu²]. Possible-pair exposures are N_rounds/2 times the product of settings. Both-nonzero X settings additionally include the published phase matching1/8 (M32/D16); vacuum settings do not. Zeng Supplement53/56/69; Zhang AppendixC basis/sifting rules.

Unreported Zhang vacuum X error counters use explicit e0=1/2 (Zeng Supplement85). Missing Zhu X00 uses the published zero-basis counter and that random-error assumption. These are conditional source-model substitutes, not measured counts. Their indices are recorded.

**Count-conditioned scope:** r_s is recovered from the observed signal count so r_p*r_s=M_Z_signal/N_rounds. Therefore pairing rate cancels in these counter-conditioned key checks; they do not independently validate Eq.4. Its closed limits are separately unit-tested and its predictions are separately compared with Zhu measured pairing rates below.

## Count-conditioned asymptotic analytical-decoy comparisons

| Source / length [km] | Calculated [bit/s] | Published finite [bit/s] | Ratio | Relative difference [%] |
| --- | --- | --- | --- | --- |
| zhang / 101 | 134781.948 | 97500 | 1.38237895 | 38.237895 |
| zhang / 202 | 9513.05687 | 6300 | 1.51000903 | 51.000903 |
| zhang / 403 | 91.0043081 | 47.8 | 1.90385582 | 90.385582 |
| zhu / 101 | 21660.9438 | 17200 | 1.2593572 | 25.935720 |
| zhu / 202 | 1569.77617 | 2080 | 0.754700081 | -24.529992 |
| zhu / 304 | 0 | 19.2 | 0 | -100.000000 |
| zhu / 407 | 0 | 0.769 | 0 | -100.000000 |

Zhang is systematically above the finite-key rates at all three lengths by 38.24–90.39%. This is a comparison of analytical asymptotic bounds with published finite-key results, not a fitted correction. Zhu has no common bias direction: +25.94%, -24.53%, then no positive analytical bound. No multiplier is introduced.

Zhu 304 km: raw C5 Y_X11 lower=-1.56868483369e-5. Zhu 407 km: raw C6/C5 error upper=0.555957968911. These do not yield positive privacy; the corresponding keys are zero versus 19.2/0.769 bit/s. A finite-key LP result and this looser analytical method are different estimators; asymptotic alone does not order their rates.

## Published finite-component arithmetic replay (Zhu only)

This check inserts the published finite M11 lower and phase-error upper into Zhu Eq.1 and TableVII signal error counts. It is **not an asymptotic prediction** or an independent decoy validation.

| Length [km] | Replay [bit/s] | Published [bit/s] | Ratio | Eq.4 pairing / measured pairing |
| --- | --- | --- | --- | --- |
| 101 | 17294.2889 | 17200 | 1.00548191 | 1.2633873 |
| 202 | 1826.79691 | 2080 | 0.878267744 | 1.06351336 |
| 304 | 19.3002481 | 19.2 | 1.00522125 | 1.02378378 |
| 407 | 0.77273735 | 0.769 | 1.00486001 | 1.03201249 |

At202 km the literal published finite components produce1826.80 rather than 2080 bit/s (-12.17%). The other arithmetic replays differ by+0.486–0.548%. No finite bound or count is altered to remove the202 km discrepancy. Zhu measured pairing uses minimum gap 63 and reference frames absent from base Eq.4, so this is not a strict identical-protocol pairing comparison.

## Zhang phase mode(b), conditional on detection gains

| Length [km] | Unfitted e_d=0 [bit/s] | Ratio | Predicted raw X | Observed weak raw X |
| --- | --- | --- | --- | --- |
| 101 | 181461.363 | 1.86114218 | 0.251625138 | 0.274267934 |
| 202 | 15072.6038 | 2.3924768 | 0.251845444 | 0.286001313 |
| 403 | 152.018565 | 3.18030472 | 0.269631178 | 0.303781898 |

Training at 101 km; only e_d fitted = 0.0928272811745.

| Length [km] | Role | Calculated [bit/s] | Ratio | Factor1.5 |
| --- | --- | --- | --- | --- |
| 101 | training | 97500 | 1 | True |
| 202 | held out | 8052.04142 | 1.27810181 | True |
| 403 | held out | 79.0997446 | 1.65480637 | False |

Training at 202 km; only e_d fitted = 0.126160038681.

| Length [km] | Role | Calculated [bit/s] | Ratio | Factor1.5 |
| --- | --- | --- | --- | --- |
| 101 | held out | 76035.4386 | 0.779850652 | True |
| 202 | training | 6300 | 1 | True |
| 403 | held out | 60.4556258 | 1.26476205 | True |

The 101 km calibration fails the inherited factor 1.5 check at 403 km (1.65480637). The 202 km calibration passes at 101/403 km (0.77985065/1.26476205). Detection gains at each held-out length are still observed inputs, not predictions. Fitted errors 9.28% and12.62% are effective discrepancy parameters, not measurements of intrinsic hardware error. They absorb finite-key/estimator, phase compensation and omitted timing differences. No hardware-only held-out validation is claimed.

## Literal source anomalies

The visually checked Zhang TableIII at 101 km prints X2mu,0=145543908 versus X0,2mu=13488265 (ratio10.7904098859). The printed value is retained in every computation, including calibration. It may affect the decoy estimate, but a typo is not assumed or repaired. Zhu TableVII labels its third column303 km; TablesI/VI and main text use304 km, which is the comparison label.

## Tests and runtime

| Check | Executed / skipped | Passed |
| --- | --- | --- |
| new_mp_tests | 13 / 0 | True |
| quick | 115 / 6 | True |
| full_units | 121 / 0 | True |

| Execution | Seconds | Exit code |
| --- | --- | --- |
| quick | 5.970072 | 0 |
| full_units | 15.370082 | 0 |
| full_T1_T7_Figure3 | 741.606271 | 0 |
| new_mp_tests | 0.411619 | 0 |
| mp_cli | 0.399802 | 0 |
| source_comparison | 0.672126 | 0 |

Quick/full unit timings above were measured during the reference T5 scan; they are not isolated benchmarks. The final MP CLI/source checks were run after that scan completed. Execution contexts are retained in details.json.

TF regression: fast 4097 / reference 65537; maximum variance error 5.10027313487e-05 (tol 1e-4), tau error 1.48396666088e-05 (tol 1e-4), key error 1.15959765034e-05 (tol 1e-3).

Full T1–T7/Figure3 reference suite passed=True; T4 retains its accepted measurement/model discrepancy. New MP tests check exact pairing limits, tiny p stability, known photon-mixture decoy inference, normalization invariance, absence of25% subtraction, worst-case privacy at 1/2, discrete averaging, spectral versus white-linewidth mode, missing inputs, units and report notices. Unit and reference logs/JSON are saved beside this report.

## Decisions and remaining limits

- Observed finite counts treated as asymptotic expectation proxies, without statistical security intervals.
- Missing vacuum X errors use explicit random e0=1/2 from Zeng Supplement85; missing Zhu X00 reuses basis0 counts.
- Literal source-table anomaly retained; no guessed typo correction.
- Zeng IID active-round gap law omits experimental minimum gaps/frame gaps; no tracking response invented.
- One source D2 error calibrated per training length; observed gains at held-out lengths remain inputs.
- Factor1.5 comparison gate transferred from the prior user-requested TF held-out check, not a physical coefficient.

Published drift inputs concern laboratory spools and do not transfer to field routes. The earlier alpha fit excludes0.5 only under its conditional residual-based95% interval; published measurement uncertainties are absent. It does not identify a correlation fraction. This warning is included in every new MP key report and USER_GUIDE.

Unverified: hardware-only detection-gain forecasts; Zhu phase forecast; a tighter source LP estimator; finite-session security intervals; real tracking residuals and experimental frame/min-gap pairing; field-route coefficients. No private or unpublished apparatus data are requested and no absent input is invented.

## Reproduce

```bash
python -m scripts.validate.validate_mp_key
python -m tfqkd.mp_key examples/mode_pairing/key_zhang202_given_gains.toml --output results/mp_key/example_zhang202
python -m scripts.report.report_mp_key
```

Primary PDFs: [Zeng](https://doi.org/10.1038/s41467-022-31534-7), Zeng supplement (`../../sources/papers/zeng2022_mode_pairing_supplement.pdf`; research archive not included in the public snapshot), [Zhang](https://doi.org/10.1103/PhysRevX.15.021037), [Zhu](https://arxiv.org/abs/2208.05649).
