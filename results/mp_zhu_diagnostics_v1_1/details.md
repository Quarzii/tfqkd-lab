# Zhu mode-pairing diagnostics — v1.1

Diagnostic only. Accepted tfqkd code, earlier tests/configuration and prior MP results are unchanged. Numbers come from executed arithmetic in audit_mp_zhu_terms.py. No new decoy method or protocol variant is enabled.

## 1. Which comparison had the −12.17% discrepancy?

It was the replay of the **published finite components**, not our analytical asymptotic decoy result. At 202 km the former is 1826.79690658 bps; the latter is 1569.77616921 bps (−24.53% versus 2080 bps). These must not be conflated.

Primary source: [Zhu PDF](https://arxiv.org/abs/2208.05649), main Eq.1/p2, denominator onp2–3; A1/p6; TableVI/p16; TableVII/p17. Rendered TableVI (`zhu_table_vi.png`; research archive not included in the public snapshot) was visually checked; M11=5.03e7 and e11=0.2571 are printed as copied.

The comparison refers to the local arXiv2208.05649v2 author PDF. [The primary version page](https://arxiv.org/abs/2208.05649) confirms v2 as the latest author version. A separate publisher supplemental file was not accessible through the web tool; its identity with the author PDF is not asserted. Version/source checks are in source_version_check.json.

### Term-by-term,202 km

| Term | Analytical/asymptotic | Finite-component replay | Published / derived | Location / interpretation |
| --- | --- | --- | --- | --- |
| r_p | 0.000456247232504 | 0.000456247232504 | 0.000429 | Zeng4 vs Zhu TableVI; replay prefactor cancels with r_s |
| r_p with source minimum gap | 0.000427907274899 | diagnostic only | 0.000429 | Zhu B7/p11; no core change |
| r_s signal fraction | 0.113536969443 | 0.113536969443 | 0.120748084138 | Not directly published; inferred from M_signal/(N*r_p) |
| r_p*r_s | 5.18009280952e-05 | 5.18009280952e-05 | 5.18009280952e-05 | Same observed TableVII signal population / N |
| M_signal | 108781949 | 108781949 | 108781949 | TableVII ZA_muB_mu total |
| q_(1,1) | 0.483058572333 | 0.462392892041 | 0.462392892041 | Published finite fraction derived from M11/M_signal; Zeng Supplement68 |
| M11 | 52548052.9796 | 50300000 | 50300000 | Our bound vs TableVI finite lower bound |
| e_(1,1)^X / phase upper | 0.278231442553 | 0.2571 | 0.2571 | Asymptotic X upper vs published finite Z-phase upper |
| e_Z | 0.000188910018518 | 0.000188910018518 | 0.000189 | Counts20550/108781949 vs rounded TableVI QBER |
| f_EC | 1.1 | 1.1 | 1.1 | TableI/p4, chosen setting 1.1 |
| N_rounds | 2.1e+12 | 2.1e+12 | 2.1e+12 | TableVI/p16 |
| N_potential_pairs | 1.05e+12 | 1.05e+12 | 1.05e+12 | Eq.1 denominator: N_rounds/2 |
| quantum rounds/s | 444830000 | 444830000 | 444830000 | TableVI header; strong/recovery already excluded |
| privacy bits | 7722996.75053 | 8936367.1358 | not directly tabulated | M11*[1-h(e11)] in Eq.1/A1 |
| error-correction bits | 312233.188421 | 312233.188421 | not directly tabulated | f*M_signal*h(e_Z), same counters and setting |
| separate additive finite penalty | 0 | 0 | none in displayed Eq.1/A1 | Finite reductions already inside M11 lower/e11 upper; intermediate penalties unreported |
| signed final key bits | 7410763.56211 | 8624133.94738 | 9819481.59971 | Published column here derived from bit/s, not an explicit counter |
| bits/potential pair | 7.05787005915e-06 | 8.21346090226e-06 | 9.34e-06 | TableVI key-rate denominator |
| bits/s | 1569.77616921 | 1826.79690658 | 2080 | Asymptotic vs finite-component replay vs TableVI final key |

### Locating the mismatch without fitting a replacement

Eq.1/A1 gives privacy=8936367.1358 bits and EC leakage=312233.188421 bits, hence K=8624133.94738 bits. The stated 2080 bps instead requires K=9819481.59971 bits; the deficit is 1195347.65233 bits. At the stated normalization and EC expense, the required privacy term is 10131714.7881 bits, 13.3762146762% above the privacy term from TableVI.

Even zero EC expense gives only 1892.93533001 bps. EC alone would require f_EC=-3.11121926282, outside f_EC>=1. An omitted nonnegative finite penalty cannot increase the replay.

The nine TableVI Sent entries sum to 2095149300000 rounds, close to the rounded2.10e12. Using that exact sent total yields 1831.02631579 bps; the discrepancy remains. The printed 9.34e-6bit/potential-pair rate correctly converts to2077.3561 bps, consistent with 2080 after printing precision. Thus the factor2 is correctly handled.

A purely algebraic one-term closure would need M11=57028235.9821 (13.3762146762% larger), or e11=0.242135935092 instead of 0.2571. Alternatively, changing only N would require 1.84436226145e+12 rounds or changing only the quantum clock would require 506485639.793 Hz. **None is a measured replacement or used by the model.** The publication does not identify which entry/provenance is responsible.

Half-last-printed-unit arithmetic envelope: replay in[1819.74740294, 1833.88194706]bps. This is an explicitly adopted deterministic printing diagnostic, not measurement uncertainty or a confidence interval. Printed precision cannot close the discrepancy.

### The actual pairing variant

Zhu B7, PDFp11, includes minimum gap 63 (Algorithm1/p7), whereas our base Zeng4 uses only maximum gap. B7 reproduces tabulated pairing rates to the rounding-level residuals below. In this count-conditioned check, r_s is inferred consistently and r_p*r_s=M_signal/N; changing the pairing prefactor therefore cannot change the replay. The final key expression A1 is the same two-term privacy-minus-EC expression, not an async vacuum-key variant.

| Length[km] | Zeng4 r_p | Zhu B7 diagnostic | TableVI r_p | B7/table−1[%] |
| --- | --- | --- | --- | --- |
| 101 | 0.00310793275591 | 0.0024598719387 | 0.00246 | -0.00520574377885 |
| 202 | 0.000456247232504 | 0.000427907274899 | 0.000429 | -0.254714475852 |
| 304 | 4.52512431572e-05 | 4.39660624951e-05 | 4.42e-05 | -0.529270373037 |
| 407 | 4.34477259931e-06 | 4.21537096611e-06 | 4.21e-06 | 0.127576392138 |

## 2. Analytical decoy versus published finite bounds, all lengths

The published q below is derived from TableVI M11 lower / TableVII signal total (Zeng Supplement68). Published e is the finite-key Z-phase upper bound, whereas ours is an asymptotic X error upper; their direct comparison diagnoses estimator tightness but is not a measurement of a pure finite-key penalty.

| Length[km] | Our q lower | Published q derived | Our q change[%] | Our phase upper | Published phase upper | Upper change[pp] | Our privacy error |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 101 | 0.519795608804 | 0.515937233641 | 0.747838091648 | 0.221223097059 | 0.2474 | -2.61769029414 | 0.221223097059 |
| 202 | 0.483058572333 | 0.462392892041 | 4.46929021786 | 0.278231442553 | 0.2571 | 2.11314425528 | 0.278231442553 |
| 304 | 0.316900179631 | 0.307692811505 | 2.99238974135 | 1 | 0.337 | 66.3 | 0.5 |
| 407 | 0.388708024544 | 0.373910465586 | 3.95751398263 | 0.555957968911 | 0.3468 | 20.9157968911 | 0.5 |

Our q lower bounds are all larger than the published finite fractions, not the cause of the zero keys. At 101 km our error upper is also tighter. At 202 km it is 2.113 pp looser; at 407 km it is 20.916 pp looser and contains the worst entropy at 1/2. At 304 km C5 is negative, so **no positive X single-photon yield is established**; the error upper is the trivial1 (not a fitted100% error), and the privacy entropy uses1/2. The published bound0.337 remains nontrivial. Exact raw yield/error bounds are preserved in JSON.

Holding our q and replacing only the phase upper by the published finite upper restores positive keys at 304/407 km. These diagnostic replacements locate the limiting estimate and are not changes to the model.

| Length[km] | Our key[bps] | Only q replaced[bps] | Only phase replaced[bps] | Both published finite bounds[bps] |
| --- | --- | --- | --- | --- |
| 101 | 21660.9438388 | 21494.1293205 | 17429.6954156 | 17294.2889037 |
| 202 | 1569.77616921 | 1499.79027124 | 1911.39768011 | 1826.79690658 |
| 304 | 0 | 0 | 20.0407578142 | 19.3002480888 |
| 407 | 0 | 0 | 0.867954159 | 0.772737349798 |

## What method does Zhu document?

- Zhang2025 C4-C7, asymptotic analytical bounds; PDF p10; plug-in observed gains
- Zhu2023 Appendix A step9/p6: decoy-state estimation, weak-weak X pairs only, with Chernoff-Hoeffding; A2/p7 and A3/p8. Ref22 points to Zeng2022.
- Zeng2022 Supplement Note3 Eqs64/65, PDF p33: linear programs; Eq67/68 p34: finite key and phase-error random sampling. The same supplement p34 notes that the programs may also be solved analytically.
- Zhu does not specify whether it used the LP or an analytical solution of the cited decoy problem, nor the intermediate single-photon X/error counts. It is not justified to label its actual solver as LP.
- Zhu EqB7/p11 and Algorithm1/p7 include L_min=63; base Zeng Eq4 does not. This explains pairing-prefactor differences but cannot alter a count-conditioned key replay.
- Zhu Eq1/p2 and A1/p6 have only privacy and f*M*h(E) terms. Chernoff/random-sampling reductions are inside tabulated M11 lower/e11 upper. No additional numerical penalty is supplied. Extra nonnegative costs cannot explain a published rate larger than the replay.

Do not assert that Zhu actually used LP: the supplied paper does not specify its numerical solver. It refers to the Zeng decoy procedure, which permits LP and an analytical solution. Zhu also states it uses only weak-weak X pairs for parameter estimation. The intermediate pre-bound X single-photon/error counts and solver choices are not supplied, so a term-by-term finite-decoy reconstruction is not identified. No method is changed to force agreement.

## Checks and decisions

Arithmetic reproduction of earlier key rates and all assertions passed; 59 model/test/config/web files retain SHA-256. Existing tests/regression are recorded in verification.json.

- Diagnostic B7 and one-term algebraic inversions are source-equation comparisons only; no new core feature.
- Half last printed unit is an engineering printing convention, not a statistical uncertainty.
- Existing count conditioning, random vacuum-error substitutions and finite/asymptotic scope retain their earlier limits.
- New clarified limit: the 202 km finite table is internally unclosed at the specified normalization; no unique erroneous scalar can be selected.
- The solver actually used by Zhu is unidentified; an LP implementation is not claimed.

## Reproduce

```bash
python -m scripts.validate.audit_mp_zhu_terms
```

## Release verification

Quick/full unit counts: quick: 115 passed, 6 skipped, full_units: 121 passed, 0 skipped.

TF fast/reference regression passed; maximum errors: {"variance_max_relative": 5.1002731348726016e-05, "tau_relative": 1.4839666608756907e-05, "key_max_relative": 1.1595976503353533e-05}.

Previous full T1–T7/Figure3 reference evidence is retained; the source scan was not repeated for documentation/diagnostics.

v1.1 updates README, STATUS, USER_GUIDE, UNKNOWNS and the supervisor summary; model, estimator and prior numerical results are unchanged.
