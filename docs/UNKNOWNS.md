# Assumptions, tool limitations, and validation limits

## Research-stage assumptions

These assumptions belong to the research reproductions in stages 1–2 and block 5.
They **do not apply to user runs**. User runs use explicit geometry, laser model,
protocol and actuator inputs, with the documented reference-index exception.

- Research comparison used two equal arms and independent stabilized lasers
  (Eq.7, bertaina2024), with published Bertaina noise and protocol coefficients.
- The research actuator scan used omega_a = 2*pi*100 .. 2*pi*100000 rad/s,
  explicitly requested as an engineering scan, not an equipment prior.
- Research scans considered fixed g and fixed g/g_crit(L,omega_a); these are
  distinct numerical studies. No such scan or length range is imposed on user inputs.
- The ideal C=g/s study has no finite first gain crossing because its delayed
  factor vanishes there. The finite actuator pole was introduced by user instruction.
- The earlier fixed tau_Q=100 ms length table was a diagnostic, not key reach;
  it remains archived in results/stage2_fixed_time/. The current core solves tau_Q.
- Historical Bertaina coefficients and source curves live in examples/ and
  the reproducibility config, not in the universal resolver's apparatus defaults.
- The measured-input experiment comparison uses user-authorized phase brackets:
  Pittaluga sigma=0..0.104 rad; Zhou sigma=0..0.099 rad, with an assumed
  nondecreasing residual RMS versus distance and, for Zhou, similar reference
  and quantum-channel RMS. These are conditional comparison hypotheses, not
  automatic assumptions for laboratory phase inputs.
- Those comparisons use effective 500 MHz quantum clocks, tau_PS=0 and published
  frame durations as working windows. This user-approved timing idealization
  avoids counting reference slots twice but does not include extra apparatus downtime.
- Out-of-sample Zhou checks fix sigma=0.05 rad and calibrate only one common e_d
  per detector projection. This fitted error also absorbs asymptotic/finite-key
  differences; it does not uniquely identify intrinsic hardware misalignment.

### Mode-pairing phase-only research extension

These entries concern `tfqkd.mp_phase` and its explicit research configuration;
they do not change TF-QKD laboratory-run assumptions or apparatus defaults.

- MP_TRACKING_RESPONSE: the approved lower band edge 1/T_track is applied only
  to divergent increment terms. It is not a tracking transfer function. An open
  publication of the estimator response/residual statistics would be needed to
  replace this band convention; no response or residual coefficient is invented.
- MP_TRANSFERRED_FIBER: l=44 and fc1=100 Hz are Bertaina Table III values from
  another line. They are used without fitting for the requested order comparison,
  not as the measured PSD of Zhou, Zhang or Zhu. Apparatus-specific open PSD data
  with normalization metadata would close this identification gap.
- MP_LASER_PSD: the curve comparison uses the explicitly stated Lorentzian
  white-frequency model (Di Domenico white case; Zhang Eq.3/Appendix B).
  A linewidth alone does not establish colored noise; an open measured laser PSD
  or an explicitly published full model would be required for an apparatus forecast.
- MP_ZHU_TRACKING_AND_BINS: Zhu Appendix B3/Fig.6 prediction is not calculated.
  Its required residual frequency-estimate statistics and within-bin distribution
  of pair intervals are not given by the published bin totals. An open residual
  characterization and interval histogram would be needed; no bin midpoint,
  distribution or estimator error is substituted.

### Mode-pairing key comparisons

These are conditional checks of `tfqkd.mp_key`, not TF-QKD user-run defaults.

- MP_OBSERVED_GAIN_CONDITIONING: Z/X gains and signal counts at every length
  are observed inputs (Zeng Supplement Eq.69). Calibrated Zhang held-out rates
  are conditional on those gains, not hardware-only out-of-sample predictions.
- MP_VACUUM_ERROR_MODEL: unreported Zhang vacuum X error counts use explicit
  e0=1/2 (Zeng Supplement Eq.85). Zhu's unreported X00 uses the published
  zero-basis counts and that same random-error assumption. Measured counts for
  these settings would close this conditional apparatus-model gap.
- MP_PAIR_INTERVAL_MODEL: Zeng Eq.4/15/16 uses IID adjacent clicks without the
  published L_min or reference-frame gaps. It is an explicit base-protocol
  approximation, not a new tracking or experimental pairing model. A complete
  published interval histogram would replace the assumed interval law.
- MP_FITTED_ERROR: one Zhang D2 error is calibrated on one length, absorbing
  finite-key/estimator and phase/timing differences; no fitted fiber coefficient,
  correction multiplier or Zhu phase forecast is introduced.

## Tool limitations

These limitations apply to user runs. No research-stage geometry or apparatus
coefficients are imported into an unspecified user installation.

| Identifier | Limitation | Interpretation / source |
|---|---|---|
| ACTUATOR_MODEL | C(s)=g/[s*(1+s/omega_a)] is the user-adopted one-pole engineering controller; Williams Eq.A6 supplies the delayed loop, not a complete real actuator. Saturation, mechanical resonances and loaded hardware dynamics are not parameterized. | Explicit pole or fitted response; no certified hardware tuning guarantee. |
| ASYMPTOTIC_KEY_MODEL | Protocols are asymptotic; finite-key effects are absent (bertaina2024 Appendix A/B). | Predicted rates do not establish a finite-session secure key. |
| MISSING_SHAPE_OR_GEOMETRY | Parametric spectral runs require explicit geometry, shape parameters and active classical pole; four-rate comparisons need both schemes. Direct nodes bypass their shape parameters. Measured-phase runs require sigma, tau, tau_PS and the selected protocol/receiver/loss data, without any spectral shapes or pole. Geometry is unnecessary there when complete arm losses are given. | Missing required inputs produce a collected list; Bertaina apparatus values are not inherited. |
| CONDITIONAL_NOISE_REQUIREMENTS | Missing r3, r2, l, s0, C2–C4 and dark counts contribute zero only for an upper key estimate. Each inverse requirement assumes the other unknown contributions are zero; individual maxima are not joint limits. | Explicit amplitudes or a separately specified joint noise budget are needed for a combined equipment guarantee. |
| N_REFERENCE_DEFAULT | n=1.4682 is a typical GROUP index at 1550 nm, Corning SMF-28 Ultra PI-1424-AEN, July 2025, p.2; not a measurement of the installed fiber. Archived Table I sensitivity over 1.44–1.48 reaches 1.818719% in variance relative to 1.45, not a universal bound. | Override the propagation index for your wavelength/fiber when known. |
| CSV_NORMALIZATION | PSD quantity, units, sidedness and pass are mandatory. Round-trip input needs an explicit calibrated PSD factor; four is not inferred for arbitrary geometry. Uncalibrated SSB dBc/Hz is rejected. | Establish normalization from estimator metadata and measurement geometry. |
| B6_MODEL_RESIDUAL | F1/Eq.6 may fail to describe the supplied spectrum. The RMS 0.5 dex acceptance gate is an engineering criterion, not a statistical agreement test. Conditional fit SE is not apparatus uncertainty. | Accepted fits remain model approximations, especially outside their measured band; no extra floor or resonance is fitted. |
| ACTUATOR_ONE_POLE | Unit-DC single-pole response fit reports conditional local SE, not calibration uncertainty. RMS magnitude defaults to 0.5 dex and RMS phase to 5 degrees; both are user-approved engineering gates, not publication values. | Thresholds are configurable; a failed or unidentifiable pole fit is rejected. |
| CLASSICAL_DETECTION_NOISE | Phase-detection noise of the classical round-trip beat measurement is not included; the supplied sources do not give its parameterization. R_classical and the classical realized ceiling fraction are optimistic estimates; comparison with R_dual is biased in favor of classical compensation. | No detection coefficient is invented. Dual retains its Appendix G detection noise; the numerical model is unchanged. |
| NONMONOTONE_INTERIOR | An imbalance range crossing zero includes that known extremum. Corners can still miss internal extrema of other nonmonotone inputs; the result is a sampled envelope. | It is not a certified continuum lower bound; another bounding method would need explicit specification. |
| MONOTONICITY_SCOPE | Table I tests and local inverse-search checks are not a proof for arbitrary protocol settings. An observed violation stops the inverse search. | Results retain their stated tested scope. |
| BAND_DIAGNOSTIC_NOT_REACH | Free fiber noise above f_b is a secondary diagnostic, not the A8 residual or a fixed-window reach criterion. Solved tau_Q can meet the phase threshold even if this ratio exceeds one. | Key benefit is calculated separately from the free-band integral. |
| SPATIAL_FIBER_MODEL | Classical remote residual assumes uniform, uncorrelated noise along each arm (Williams Eq.A8/A11). Fiber model coefficients scale linearly with arm length (bertaina2024 Eq.6/8); an uploaded spectrum does not separately identify spatial noise correlations. | Predictions use these source-model assumptions; they do not certify another noise distribution or an extrapolated route. |
| SHARED_LASER_MODEL | The existing independent-laser core uses the same PSD for both lasers; shared controller gain/pole is assumed for the two classical arms. Applying Williams remote residual to common-laser Eq.5 with the original K is a documented engineering composition. | Different laser PSDs or different arm controllers require an explicitly specified extension; do not infer their identity from a hardware name. |
| SNS_SIGNAL_COUPLING | The SNS wrapper fixes signal intensity s=decoy_big/2 and not-send intensity n=decoy_mini/2, following the Bertaina author plotting configuration. They are not independent user inputs. | Bertaina Appendix D defines decoy_big/medium/mini as TOTAL two-user intensities; per-user symmetric values must be multiplied by two. Zhou2023 symmetric Table S4 is compatible: total decoy2=0.986 and decoy0=0.0004 give s=0.493 and n=0.0002. Independent signal intensities or two asymmetric source sets remain unsupported; no tool extension was made during the audit. |
| SCALAR_RECEIVER_AND_LOSS_MODEL | The core protocol uses scalar detector efficiency/noise/error. The wrapper accepts complete arm losses including insertions, with author equalization effective loss=2*max(loss_A,loss_B). Two explicit detectors produce two scalar receiver projections; no detector averaging is performed. | The smaller projection is a displayed reference, not a certified bound for two unequal detectors. Explicit dark/background counts are retained separately and their supplied total is used in pdc (Zhou2023 Supplementary Note3/Table S3). |
| MEASURED_PHASE_SCOPE | A supplied residual sigma and operating window bypass the spectral prediction. Eq.1 small-phase error is reported; the protocol retains its unchanged author Gaussian phase-error function. The input does not identify noise components, other compensation schemes, a loop gain or attainable reach at another length. | It is a conditional key estimate from an observation or explicitly labelled phase bound, not a noise forecast. Missing intrinsic error/fEC use user-authorized ideal limits, with an upper-model-estimate label. |
| DIRECT_PSD_SCOPE | Direct PSD is inserted without a parametric fit. Log-log interpolation and explicitly selected endpoint power-law extrapolation are numerical engineering choices. Measurements outside the supplied frequency band are not created by interpolation. | Extrapolation is rejected by default and flagged when enabled. No inverse coefficient requirements or coefficient sensitivity for the direct node; mode=fit remains available. |
| DIRECT_ARM_STATE | Uploaded single-arm fiber PSD describes the selected compensation state, excluding laser/central detection contributions. It does not identify another compensation state or a changed controller. Linear arm-length scaling uses the Eq.6/8 source-model assumption, not a calibrated route map. | A direct arm run omits four-scheme ceiling and loop tuning. Use separately measured states or a fitted source model for such predictions. |
| CONDITIONAL_PROTOCOL_REQUIREMENTS | Unreported e_d/f_EC use ideal limits only for a labelled upper model estimate. An explicit target rate is needed for inverse requirements; e_d limits fix f_EC=1.16 (Liu2023 Sec.2.1, Eq.1 discussion), f_EC limits condition on the supplied or ideal e_d. | Separate conditional limits per scalar projection, not joint equipment guarantees. No target rate is invented. |
| REACH_IDENTIFIABILITY | Reach scans are disabled by default; enable reach.enabled explicitly. A phase/window observation or complete-loss point does not supply a length dependence; reach remains unidentified. Spectral/per-km runs retain fixed imbalance (the existing Part C engineering convention), using an automatic numerical bracket from the input length or explicit bounds. Target-rate shifts use actual numerical crossings, which need not be symmetric. | A multiplier does not move the mathematical zero. Finite bracketing cannot certify hidden crossings; failed automatic brackets and crossings outside explicit bounds are reported as unidentified. |


## Validation limits

These are limits of comparisons with the supplied publications and datasets.
They are recorded here and in validation reports, not reproduced as the Tool
limitations section of each user report.

| Identifier | Validation limit | Available evidence |
|---|---|---|
| T4_UNCERTAINTY | В отдельном измерительном TXT нет неопределённостей и числа усреднений. Расхождение Eq. 6 по уровню и форме ненулевое. | Сохраняются численные невязки; статистическая значимость не утверждается. Это предел валидации, а не незавершённая реализация. |
| T4_HIGH_FREQUENCY | Вклады резонансов и измерительного фона не разделены для этого TXT. | Возможные причины отмечены как гипотезы; l и fc1 не изменяются. |
| APPARATUS_IDENTITY | Разрешённый исполнительный полюс не идентифицирован с конкретной экспериментальной аппаратурой. | Результаты представлены как функции omega_a и g; точное совпадение с аппаратным servo bump Williams не заявляется. |
| B6_RMS_CLOSURE | Оцифрованный спектр Jiang даёт около 9.920 rad на доступной полосе, текст сообщает 19.2 rad для 1–1000 Hz. Нормировка установлена по Fig.1/3 и не меняется ради замыкания RMS. | Открытые исходные данные/пояснение этого расхождения. Это ограничение валидации, не запрос личных сведений. |
| B6_LASER_IDENTIFIABILITY | В архивных спектрах Snigirev fc F1 не определяется; Orbits имеет вырождение r2/fc. NF6700 и Droste отклонены по невязке. Отклонённые коэффициенты не применяются к расчёту ключа. | Измерительный диапазон, определяющий все параметры F1/Eq.6, и установленная нормировка. |
| EXPERIMENT_END_TO_END_UNVERIFIED | The spectral audit did not establish an end-to-end prediction of the measured key rates of Pittaluga2021 or Zhou2023. The new measured-phase path executes conditional comparisons using explicit phase brackets and timing idealizations, not same-window quantum-channel measurements at each length. Clivati2022 did not execute a full QKD transmission. | All four published rates lie outside the tested conditional phase/projection envelopes; parameters are not adjusted to force agreement. See MEASURED_INPUTS_REPORT.md and results/measured_inputs/experiment_comparison.json. |
| PHASE_SOURCE_ATTRIBUTION | Prior audit assignments of Zhou sigma=0.099/0.096/0.097 to three lengths are withdrawn. Figure3e–g describes frequency-offset distributions on lambda_c at615.6km. | Original results are preserved under results/audit/withdrawn_phase_attribution. Current comparisons use explicitly user-approved phase brackets, not those erroneous per-length measurements. |
| HELD_OUT_PROTOCOL_EDGE | Calibrating one common e_d per projection gives held-out ratios 0.8902–1.0786 between Zhou's two shorter lengths, within factor1.5. The calibrated 615.59km rate is zero versus 0.32bit/s published. | The scope is the unchanged asymptotic protocol/loss model plus a fitted error that absorbs finite-key differences; no general multiplicative correction is identified. |
| PITTALUGA_CAL_QBER | No CAL QBER or numeric f_EC located in the supplied CAL table/section; fitted e_d is conditional on ideal f_EC and the approved phase bracket. | CAL plausibility against its own published QBER remains unverified; SNS tables are not interchangeable. |
| MP_GRAPHIC_COMPARISON | Published vector traces/marker centers supply coordinates, not raw event covariances or a statistical agreement test. The unfitted transferred fiber PSD gives nonzero drift and X-error discrepancies. | Full numbers, digitization calibration, source-model controls and conditional assumptions are in results/mp_phase/details.md and JSON; phase validation remains distinct from the new MP key checks. |
| MP_ANALYTICAL_DECOY_LIMIT | Zhang C4-C7 is an allowed analytical method; the supplied Zhu paper does not identify its actual decoy solver. For Zhu304km, C5 is negative; at407km the error upper bound exceeds1/2. Both establish zero key with this estimator, despite Z single-photon fractions above the published finite fractions. | Zhu Appendix A documents decoy estimation plus Chernoff-Hoeffding and cites Zeng, whose Supplement64/65 allows LP and also an analytical solution. Do not identify Zhu's actual implementation as LP without evidence. No estimator is changed. |
| MP_COUNT_PLUGIN | Finite observed counters are used as asymptotic expectation proxies without confidence intervals. | These are arithmetic/model comparisons, not finite-session secure-key claims. |
| MP_LITERAL_TABLES | Zhang101km prints X2mu,0=145543908 versus X0,2mu=13488265 (10.7904 ratio). Zhu202km finite-component replay gives1826.80bps versus2080bps (-12.17%). | Values and residuals are retained literally; presumed typographical corrections are not used. |
| MP_ZHU_202_CLOSURE | The printed finite M11/e11 and stated normalization imply insufficient privacy bits at202km. Even zero EC expense reaches only1892.94bps; exact sent totals and printing precision do not close2080bps. | Zhu Eq1/A1 matches the two-term key expression. B7 corrects the minimum pairing gap but cancels with r_s under count conditioning. A unique erroneous table entry/provenance cannot be established. See results/mp_zhu_diagnostics_v1_1/. |
| MP_HELD_OUT_SCOPE | Calibration on Zhang202km gives conditional held-out ratios0.780/1.265; calibration on101km gives1.655 at403km, outside factor1.5. | Conditional on observed detection gains; no independent hardware-only forecast or universal correction factor is established. |

Historical reports and digitization remain in results/, B6_SOURCE_AUDIT.md, and examples/b6/.
