# Out-of-sample protocol check and target-rate wrapper

## Scope and unchanged calculations

Only e_d is calibrated. Spectra, integration, Gaussian author phase-error function, scalar protocols, losses, intensities, pulse clocks and apparatus coefficients are unchanged. The 16 historical kernel/config files pass the existing path-aware identity check. Reorganization was already completed and accepted; it is not repeated.

Calibration does not identify intrinsic misalignment uniquely: it also absorbs the asymptotic/finite-key discrepancy. No fitted rate multiplier is used.

## Published inputs and source locations

- Zhou2023: sources/papers/zhou2023.pdf; Tables S1/S2 fiber lengths/losses, Note3/Table S3 receiver efficiency, intrinsic dark and scattered counts and Charlie insertion losses, Table S4 symmetric per-user intensities/probabilities, Methods fEC=1.1, Table S5 measured finite-key rates and X11 QBER. Exact resolved inputs are the preserved examples/measured_inputs/zhou2023_* files. Total two-user decoys and author signal coupling follow Appendix D / QKD.ipynb in Bertaina, not imputed apparatus coefficients.
- Pittaluga2021: sources/papers/pittaluga2021.pdf; Tables I/II receiver/channel transmission, Table III CAL signal and 852.7 bit/s, Note VI code-basis normalization and pattern length. No own CAL QBER or numeric fEC located; SNS Tables V-X are not substituted. CAL numerical truncation comes from QKD.ipynb.
- Conditional phase bounds: Pittaluga Fig2f (605km), Zhou Fig3e-g (615.6km, lambda_c channel, frequency offsets). User-approved nondecreasing residual RMS with length and similar Zhou reference/quantum-channel RMS; not same-window measurements at all compared lengths.
- Effective 500MHz clock already includes reference slots. User-approved tau_PS=0; the published pattern periods are used as windows. This is a timing idealization, not a new observation.
- For Zhou fits, sigma=0.05rad is fixed by the user at all lengths. For Pittaluga, both previously approved phase endpoints are fitted separately; fEC remains the labelled ideal value 1 because its apparatus value was not published.
- Conditional e_d requirements only: fEC=1.16, Liu et al., Quantum Frontiers 2,16 (2023), Sec2.1 discussion of Eq1, https://link.springer.com/article/10.1007/s44214-023-00039-9 . Primary-source excerpt and scope saved in sources/data/f_ec_reference.txt. No claim that it is a universal standard or the apparatus value.

## Zhou: calibration and held-out predictions

| Anchor km | Projection | Fitted e_d | Test km | Published bit/s | Model bit/s | Signed bound bit/s | Ratio | Held out | Protocol edge | Within factor1.5 | Below X11 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 403.73 | D0 | 0.03522196347503552 | 403.73 | 146.7 | 146.69999999999996 | 146.69999999999996 | 0.9999999999999998 | False | False | True | True |
| 403.73 | D0 | 0.03522196347503552 | 518.16 | 14.38 | 12.853403561748456 | 12.853403561748456 | 0.8938389124998926 | True | False | True | True |
| 403.73 | D0 | 0.03522196347503552 | 615.59 | 0.32 | 0.0 | -3.347049833808404 | 0.0 | True | True | False | True |
| 403.73 | D1 | 0.03962004967709673 | 403.73 | 146.7 | 146.69999999999828 | 146.69999999999828 | 0.9999999999999883 | False | False | True | True |
| 403.73 | D1 | 0.03962004967709673 | 518.16 | 14.38 | 12.801008627894838 | 12.801008627894838 | 0.8901953148744671 | True | False | True | True |
| 403.73 | D1 | 0.03962004967709673 | 615.59 | 0.32 | 0.0 | -3.4354150661170877 | 0.0 | True | True | False | True |
| 518.16 | D0 | 0.03140869248500018 | 403.73 | 146.7 | 157.8560941755968 | 157.8560941755968 | 1.0760469950620097 | True | False | True | True |
| 518.16 | D0 | 0.03140869248500018 | 518.16 | 14.38 | 14.38 | 14.38 | 1.0 | False | False | True | True |
| 518.16 | D0 | 0.03140869248500018 | 615.59 | 0.32 | 0.0 | -3.1256389381986285 | 0.0 | True | True | False | True |
| 518.16 | D1 | 0.03571252144474864 | 403.73 | 146.7 | 158.22675866641944 | 158.22675866641944 | 1.0785736787076992 | True | False | True | True |
| 518.16 | D1 | 0.03571252144474864 | 518.16 | 14.38 | 14.380000000000008 | 14.380000000000008 | 1.0000000000000004 | False | False | True | True |
| 518.16 | D1 | 0.03571252144474864 | 615.59 | 0.32 | 0.0 | -3.2049349085140877 | 0.0 | True | True | False | True |

The common e_d for each projection satisfies the user QBER criterion at every length. The shortest-to-middle predictions underpredict by about11%; reverse predictions overpredict by about8%. These opposite signs do not identify a single calibrated systematic multiplier. All calibrated longest-point predictions are zero versus0.32bit/s published; raw negative bounds are preserved above. No correction factor is introduced.

## Pittaluga: single-point conditional calibration

| Projection | sigma rad | e_d | fEC condition | Fitted bit/s | Relative fit residual | Own CAL QBER |
| --- | --- | --- | --- | --- | --- | --- |
| D0 | 0.0 | 0.031743596545887307 | 1.0 | 852.7000000817238 | 9.584111282379126e-11 | unverified: no CAL QBER located in Table III or CAL supplementary section; SNS errors in Tables V-X are not substituted |
| D0 | 0.104 | 0.029199793614830053 | 1.0 | 852.6999998951678 | 1.2294154583258887e-10 | unverified: no CAL QBER located in Table III or CAL supplementary section; SNS errors in Tables V-X are not substituted |
| D1 | 0.0 | 0.03593895466119464 | 1.0 | 852.6999999004444 | 1.1675349576023564e-10 | unverified: no CAL QBER located in Table III or CAL supplementary section; SNS errors in Tables V-X are not substituted |
| D1 | 0.104 | 0.03341794221624105 | 1.0 | 852.6999999873855 | 1.479372180313021e-11 | unverified: no CAL QBER located in Table III or CAL supplementary section; SNS errors in Tables V-X are not substituted |

There is no held-out CAL point. A root matching the one published rate is not validation. The requested e_d/QBER plausibility comparison remains unresolved because no CAL QBER was located; it is not replaced by an SNS error or a phase-only contribution.

## Phase-bound endpoint effects

Percent changes are recomputed at a fixed fitted e_d; they are not a new fit at each phase endpoint. A zero baseline gives an undefined percentage.

| Paper/profile | Length km | Projection | e_d | R(sigma=0) | R(sigma=upper) | Change % |
| --- | --- | --- | --- | --- | --- | --- |
| Zhou, anchor 403.73 | 403.73 | D0 | 0.03522196347503552 | 148.47872296260903 | 141.62132459390702 | -4.618438407790515 |
| Zhou, anchor 403.73 | 518.16 | D0 | 0.03522196347503552 | 13.096846208851993 | 12.158233845752054 | -7.166705236758009 |
| Zhou, anchor 403.73 | 615.59 | D0 | 0.03522196347503552 | 0.0 | 0.0 | None |
| Zhou, anchor 403.73 | 403.73 | D1 | 0.03962004967709673 | 148.49414042147782 | 141.5733792540601 | -4.660629131745009 |
| Zhou, anchor 403.73 | 518.16 | D1 | 0.03962004967709673 | 13.04681146276918 | 12.09858300348925 | -7.2678942436305345 |
| Zhou, anchor 403.73 | 615.59 | D1 | 0.03962004967709673 | 0.0 | 0.0 | None |
| Zhou, anchor 518.16 | 403.73 | D0 | 0.03140869248500018 | 159.7536203977304 | 152.4428326338368 | -4.576289254473432 |
| Zhou, anchor 518.16 | 518.16 | D0 | 0.03140869248500018 | 14.639590197824733 | 13.639332033539818 | -6.832555766714988 |
| Zhou, anchor 518.16 | 615.59 | D0 | 0.03140869248500018 | 0.0 | 0.0 | None |
| Zhou, anchor 518.16 | 403.73 | D1 | 0.03571252144474864 | 160.13832465751528 | 152.76826937147626 | -4.6023057265031415 |
| Zhou, anchor 518.16 | 518.16 | D1 | 0.03571252144474864 | 14.641805212743 | 13.63233278923494 | -6.894453305726945 |
| Zhou, anchor 518.16 | 615.59 | D1 | 0.03571252144474864 | 0.0 | 0.0 | None |
| Pittaluga, fit sigma 0.0 | 368.702 | D0 | 0.031743596545887307 | 852.7000000817238 | 804.4498942133608 | -5.658508955522301 |
| Pittaluga, fit sigma 0.104 | 368.702 | D0 | 0.029199793614830053 | 902.3170712128336 | 852.6999998951678 | -5.498851002671811 |
| Pittaluga, fit sigma 0.0 | 368.702 | D1 | 0.03593895466119464 | 852.6999999004444 | 803.9328472800922 | -5.71914537657393 |
| Pittaluga, fit sigma 0.104 | 368.702 | D1 | 0.03341794221624105 | 902.7635991288927 | 852.6999999873855 | -5.545593463207343 |
| zhou2023, ideal e_d | 403.73 | D0 | 0 | 304.97752296077596 | 284.69389680070503 | -6.650859369291839 |
| zhou2023, ideal e_d | 403.73 | D1 | 0 | 331.25921599992995 | 309.2498809394249 | -6.64414271285051 |
| zhou2023, ideal e_d | 518.16 | D0 | 0 | 34.17116565240214 | 31.606518117422436 | -7.505297188477433 |
| zhou2023, ideal e_d | 518.16 | D1 | 0 | 37.72433711230537 | 34.92892816583992 | -7.410094279837221 |
| zhou2023, ideal e_d | 615.59 | D0 | 0 | 0.0 | 0.0 | None |
| zhou2023, ideal e_d | 615.59 | D1 | 0 | 0.03321076615838136 | 0.0 | -100.0 |
| pittaluga2021, ideal e_d | 368.702 | D0 | 0 | 1605.5440723524218 | 1521.530140335468 | -5.232739073543935 |
| pittaluga2021, ideal e_d | 368.702 | D1 | 0 | 1733.264499799295 | 1644.0009428041062 | -5.150025111892909 |

For positive calibrated rates, the endpoint change is4.58–7.27% (Zhou) and5.50–5.72% (Pittaluga). These comparisons test primarily the protocol/loss bookkeeping under conditional phase assumptions, not a spectral phase-noise forecast. At the uncalibrated Zhou615.59km D1 edge the small positive0.03321bit/s goes tozero (100%); D0 and both calibrated projections have undefined relative changes. The bulk percentage is not claimed for that edge.

## Wrapper requirements and reach demonstration

User runs with missing e_d/fEC now start with an upper-model-estimate label in Markdown and HTML. No apparatus value or target rate is silently supplied. requirements.target_key_bps explicitly sets the inverse-search goal; absent targets give an explanatory status. Separate detector projections receive separate conditional requirements, not a joint guarantee.

| Parameter | Projection | Maximum | Conditional fEC | Conditional e_d | Verified bit/s | Target relative residual | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| keyrate.detector_error | 0 | 0.03456567698967732 | 1.16 | None | 146.7000065101102 | 4.4377029384534694e-08 | conditional single-parameter requirement |
| keyrate.detector_error | 1 | 0.03893022556563227 | 1.16 | None | 146.7000448482087 | 3.0571376075094747e-07 | conditional single-parameter requirement |
| keyrate.f_error | 0 | 5.981454576626453 | None | 0.0 | 146.70000000000002 | 2.220446049250313e-16 | conditional single-parameter requirement |
| keyrate.f_error | 1 | 6.4822815102381 | None | 0.0 | 146.69999999999996 | 2.220446049250313e-16 | conditional single-parameter requirement |

Demonstration inputs are saved as protocol_target_input.json and reach_input.json. Reports and all provenance are in protocol_target_example/ and reach_example/. The spectral reach example explicitly selects the Bertaina fixture, fixed imbalance, and a target equal to its own calculated rate at400km; it is not a forecast for an identified laboratory.

Reach JSON: {"zero_rate_reach_km": 547.547229349694, "target_key_bps": 823.3471421044576, "input_case_calculations": 52, "spectral_calculations": 52, "message": "A factor of 2 in key rate corresponds to 27.9155787 km (doubling) / 28.8325273 km (halving) in reach for this configuration.", "rate_curve_basis": "explicit point inputs and their labelled ideal limits, if any", "status": "evaluated on source-model length scan", "geometry": "fixed_imbalance", "domain_total_km": [100.0, 1000.0], "domain_source": "explicit user numerical bounds", "zero_rate": {"threshold_bps": 0.0, "roots_km": [547.547229349694], "reach_km": 547.547229349694, "status": "crossing found within explicit length domain", "multiple_crossings": false}, "sampled_curve": {"length_total_km": [100.0, 128.125, 156.25, 184.375, 212.5, 240.625, 268.75, 296.875, 325.0, 353.125, 381.25, 409.375, 437.5, 465.625, 493.75, 521.875, 550.0, 578.125, 606.25, 634.375, 662.5, 690.625, 718.75, 746.875, 775.0, 803.125, 831.25, 859.375, 887.5, 915.625, 943.75, 971.875, 1000.0], "key_bps": [906531.31416707, 468658.8867882344, 243533.358721022, 126885.31705155605, 66191.28571490849, 34541.19534997407, 18017.827806959158, 9387.014864095683, 4878.258354084498, 2523.491204357302, 1294.6324050136575, 654.3923078450623, 321.9047653149839, 150.3252344399862, 62.86786987363373, 19.367601124868038, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "raw_key_bps": [906531.31416707, 468658.8867882344, 243533.358721022, 126885.31705155605, 66191.28571490849, 34541.19534997407, 18017.827806959158, 9387.014864095683, 4878.258354084498, 2523.491204357302, 1294.6324050136575, 654.3923078450623, 321.9047653149839, 150.3252344399862, 62.86786987363373, 19.367601124868038, -1.2078278477737576, -9.910335278030262, -12.60274087320615, -12.452342352361752, -8.922576893901514, -6.084611049429871, -4.292136181379801, -3.2789148919270668, -2.7347869785609586, -2.448501307793518, -2.29881815167523, -2.22061715381264, -2.179738831667398, -2.1583594779011115, -2.1471748862329236, -2.141322836985513, -2.138260685550888]}, "limitation": "Scan brackets sampled crossings; hidden internal crossings are not certified. Same apparatus coefficients and imbalance; automatic g is re-optimized at each length. Uploaded arm length scaling uses the documented Eq.6/8 model assumption.", "target_reach": {"nominal": {"threshold_bps": 823.3471421044576, "roots_km": [400.0000009283513], "reach_km": 400.0000009283513, "status": "crossing found within explicit length domain", "multiple_crossings": false}, "doubled_rate": {"threshold_bps": 411.6735710522288, "roots_km": [427.915579597391], "reach_km": 427.915579597391, "status": "crossing found within explicit length domain", "multiple_crossings": false}, "halved_rate": {"threshold_bps": 1646.6942842089152, "roots_km": [371.1674736115781], "reach_km": 371.1674736115781, "status": "crossing found within explicit length domain", "multiple_crossings": false}}, "doubling_shift_km": 27.91557866903969, "halving_shift_km": 28.832527316773223}

Automatic-bracket verification: {"zero_rate_reach_km": 547.5473437156237, "target_key_bps": 823.3471421044576, "input_case_calculations": 51, "spectral_calculations": 51, "message": "A factor of 2 in key rate corresponds to 27.9156497 km (doubling) / 28.8324532 km (halving) in reach for this configuration.", "rate_curve_basis": "explicit point inputs and their labelled ideal limits, if any", "status": "evaluated on source-model length scan", "geometry": "fixed_imbalance", "domain_total_km": [200.0, 799.94], "domain_source": "automatic numerical bracket from input length; same fixed-imbalance convention as Part C", "zero_rate": {"threshold_bps": 0.0, "roots_km": [547.5473437156237], "reach_km": 547.5473437156237, "status": "crossing found within explicit length domain", "multiple_crossings": false}, "sampled_curve": {"length_total_km": [200.0, 218.74812500000002, 237.49625, 256.244375, 274.9925, 293.740625, 312.48875, 331.236875, 349.985, 368.73312500000003, 387.48125000000005, 406.229375, 424.9775, 443.72562500000004, 462.47375, 481.221875, 499.97, 518.7181250000001, 537.4662500000001, 556.214375, 574.9625000000001, 593.710625, 612.45875, 631.2068750000001, 649.955, 668.703125, 687.4512500000001, 706.199375, 724.9475, 743.6956250000001, 762.44375, 781.1918750000001, 799.94], "key_bps": [88382.71031835252, 57286.0291582563, 37133.20951772852, 24066.708983947097, 15592.374146209113, 10095.486739680686, 6529.785602583184, 4216.959405438301, 2717.078670372242, 1744.7466574767056, 1114.7845743994988, 707.0253111563871, 443.48296410324485, 273.5432185540336, 164.35458483997235, 94.59301416007318, 50.414971192303966, 22.83037308138822, 5.996742520956475, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "raw_key_bps": [88382.71031835252, 57286.0291582563, 37133.20951772852, 24066.708983947097, 15592.374146209113, 10095.486739680686, 6529.785602583184, 4216.959405438301, 2717.078670372242, 1744.7466574767056, 1114.7845743994988, 707.0253111563871, 443.48296410324485, 273.5432185540336, 164.35458483997235, 94.59301416007318, 50.414971192303966, 22.83037308138822, 5.996742520956475, -3.888920325977432, -9.310215430695994, -11.898861142214722, -12.736717225500716, -12.548121628853586, -10.492681527792376, -8.2057829068178, -6.3502093112337015, -4.975356394642497, -4.013916064104859, -3.3649325584943286, -2.935730571390904, -2.655013695262145, -2.4724140365398473]}, "limitation": "Scan brackets sampled crossings; hidden internal crossings are not certified. Same apparatus coefficients and imbalance; automatic g is re-optimized at each length. Uploaded arm length scaling uses the documented Eq.6/8 model assumption.", "target_reach": {"nominal": {"threshold_bps": 823.3471421044576, "roots_km": [399.99992286385077], "reach_km": 399.99992286385077, "status": "crossing found within explicit length domain", "multiple_crossings": false}, "doubled_rate": {"threshold_bps": 411.6735710522288, "roots_km": [427.91557258492605], "reach_km": 427.91557258492605, "status": "crossing found within explicit length domain", "multiple_crossings": false}, "halved_rate": {"threshold_bps": 1646.6942842089152, "roots_km": [371.16746969534444], "reach_km": 371.16746969534444, "status": "crossing found within explicit length domain", "multiple_crossings": false}}, "doubling_shift_km": 27.915649721075283, "halving_shift_km": 28.832453168506333}

Measured-phase example: not identified: one measured phase/window does not supply phase versus length. No length or loss extrapolation is inserted. A rate multiplier cannot shift the mathematical zero; the user-approved kilometre shifts refer to target-rate crossings. The shifts are asymmetric and calculated through the actual kernel rather than an analytic loss law.

Timing: phase with four conditional protocol requirements=0.3446256160386838s; spectral report with reach, sensitivity and four-scheme ceiling=18.113601811986882s. These are measured executions, not a guaranteed latency for all equipment/settings.

## Verification gates

| What | How | Result |
| --- | --- | --- |
| Fast unit tests | unittest discover -s tests | 73 passed; 154.288999656972s |
| Fast/reference | Table I grid4097 versus65537, unchanged calibration gate | {"points": 4097, "variance_max_relative": 5.1002731348726016e-05, "tau_max_relative": 1.4839666608756907e-05, "key_max_relative": 1.1595976503353533e-05, "passed": true} |
| Direct PSD path | 14 cases on fast/reference grids | 2.220446049250313e-16 maximum relative error |
| Parallel schedule | Existing ProcessPoolExecutor identity gate | passed; sandbox IPC denied once, same command rerun outside sandbox |
| Inverse protocol target | Returned maximum fed back into unchanged kernel | four projection/parameter combinations within1e-3 relative target tolerance |
| Distance roots | Kernel reevaluated at target/half-target/twice-target roots; explicit/automatic bracket comparison | passed |
| Physics identity | 16 existing baseline hashes, path-only relocation normalization | no changed files |
| Full reference T1–T7 | Previously executed suite; kernel identity reconfirmed | not rerun for wrapper-only additions; previous accepted T4 discrepancy retained |
| Zhou interior held out | one e_d fit, predict other shorter length; factor1.5 | passed |
| Zhou615.59km | same fitted e_d | failed: zero versus0.32bit/s |
| Pittaluga plausibility | needed e_d versus its own CAL QBER | unverified: CAL QBER not located |

## Decisions and remaining limits

- One fitted e_d per detector projection absorbs finite-key/model differences; fixed sigma=0.05 is user-authorized.
- Pittaluga calibrated at both approved phase endpoints with ideal unpublished fEC; not an out-of-sample validation.
- Operational factor-two distance uses the user-approved target/half-target/twice-target inversions; no analytic attenuation approximation.
- Fixed-imbalance length scan reuses the existing Part C engineering convention; automatic numerical brackets start from the input length.
- Conditional e_d requirement fixes literature example fEC=1.16; main calculation retains explicit or labelled ideal fEC.
- No multiplicative correction is introduced: the longest Zhou point fails and interior calibrated ratios straddle unity.

No changes to classical detection-noise omission, scalar detector projection, asymptotic key, T4 discrepancy or measurement normalization. Tool limitations add conditional protocol requirements and reach-identifiability rules; the current user report includes only Tool limitations. The user guide documents the tested systematic deviations and their restricted scope. Part D remains deferred.
