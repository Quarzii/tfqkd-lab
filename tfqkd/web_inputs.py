"""Web-only input descriptions and preflight validation; no model defaults.

Bounds and unusual-value thresholds are user-requested interface policy.
SNS intensity normalization: Appendix D, bertaina2024 / QKD.ipynb Cell 21.
"""
import math

# Display schema shared by the form and the server's preflight checks.
FIELDS = {}

def field(id, path, label, help, rule=None, unusual=None):
    FIELDS[id] = dict(path=path, label=label, help=help, rule=rule, unusual=unusual,
                      variable=path or id)

field('sigma','phase.sigma_phi_rad','Residual phase noise RMS [rad], measured over window tau','Observed residual RMS on the specified window; this mode estimates a measured operating point, not a spectral forecast.','phase','high_phase')
field('tau','phase.tau_s','Key window tau [s]','Window associated with the reported residual phase RMS.','positive')
field('taups','phase.tau_ps_s','Time spent on phase stabilization between key windows [s]','Additional downtime only; do not count reference slots twice if the clock already excludes them.','nonnegative')
for id, path, label, help in [
 ('length','line.length_km','Total physical fiber length [km]','Sum of the two arm lengths.'),
 ('lossA','line.loss_a_db','Complete arm A loss [dB]','Total measured loss including insertions; supply both arms.'),
 ('lossB','line.loss_b_db','Complete arm B loss [dB]','Total measured loss including insertions; overrides per-km attenuation when both arms are supplied.'),
 ('alpha','line.attenuation_db_per_km','Fiber attenuation [dB/km]','Used only when complete arm losses are not supplied.'),
 ('clockrate','keyrate.clockrate_hz','Effective quantum pulse rate [Hz]','Pulse rate used by the protocol; exclude reference slots only once.'),
 ('targetKey','requirements.target_key_bps','Required key rate [bit/s]','Optional target for conditional equipment requirements and target reach.')]:field(id,path,label,help,'positive')
field('imbalance','line.imbalance_km','Arm length difference A − B [km]','Signed difference; its magnitude must be less than the total length.','finite')
field('fError','keyrate.f_error','Error correction inefficiency f_EC [dimensionless]','Ratio to the Shannon limit; at least 1. Blank means an explicitly labelled ideal upper estimate.','fec','high_fec')
for id, path, label in [('detectorEfficiency','detector.efficiency','Detector efficiency'),('d0Efficiency','detector.channels.0.efficiency','D0 efficiency'),('d1Efficiency','detector.channels.1.efficiency','D1 efficiency')]:
    field(id,path,label+' [fraction]','Detection probability in (0, 1]. Two detectors use separate scalar model projections, not an average.','probability','low_efficiency')
for id,path,label in [('detectorDark','detector.dark_count_rate_hz','Detector dark counts'),('d0Dark','detector.channels.0.dark_count_rate_hz','D0 dark counts'),('d1Dark','detector.channels.1.dark_count_rate_hz','D1 dark counts'),('d0Background','detector.channels.0.background_count_rate_hz','D0 background counts'),('d1Background','detector.channels.1.background_count_rate_hz','D1 background counts')]:
    field(id,path,label+' [Hz]','Nonnegative observed noise-count rate; background is added to dark counts in the scalar projection.','nonnegative')
for id in ('detectorError','twoDetectorError'):
    field(id,'detector.error','Intrinsic misalignment error e_d [fraction]','Non-phase intrinsic error, in [0, 0.5). Blank produces an upper estimate with conditional requirements.','error')
for id,key,label in [('decoyBig','decoy_big','Signal intensity'),('decoyMedium','decoy_medium','Medium decoy intensity'),('decoyMini','decoy_mini','Weak decoy intensity')]:
    field(id,'keyrate.'+key,label+' [photons/pulse/user] — per user, as usually published','SNS: enter one-user intensity. The form multiplies by 2 for the author core convention (Bertaina Appendix D); both values are reported.','nonnegative')
field('pzSns','keyrate.pz_sns','Probability of the code (Z) basis [fraction]','SNS probability in (0, 1].','probability')
field('epsSns','keyrate.eps_sns_aopp','SNS sending probability [fraction]','SNS requires both sending and not-sending events; epsilon must be strictly between 0 and 1. Published values are typically 0.05-0.3.','sns_probability','sns_typical')
field('pzCal','keyrate.pz_cal','Probability of the code (Z) basis [fraction]','CAL probability in (0, 1].','probability')
field('uCal','keyrate.u_cal','CAL code intensity [photons/pulse/user] — per user, as usually published','CAL α² is already a one-user intensity in Bertaina Appendix B / QKD.ipynb Cell 19.','positive')
field('nminCal','keyrate.nmin_cal','CAL photon-number lower bound [integer]','Lower summation truncation, nonnegative integer.','integer')
field('nmaxCal','keyrate.nmax_cal','CAL photon-number upper bound [integer]','Upper summation truncation, integer greater than the lower bound.','integer')
for id,path,label,help,rule in [
 ('laserLinewidth','laser.lorentz_width_hz','Lorentzian linewidth [Hz]','White frequency-noise contribution only; does not identify flicker noise.','positive'),
 ('r3','laser.r3','Laser flicker coefficient [rad² Hz²]','F1 coefficient; blank means an unmeasured amplitude and an optimistic zero contribution.','nonnegative'),
 ('r2','laser.r2','Laser white-frequency coefficient [rad² Hz]','F1 coefficient for the 1/f² phase term.','nonnegative'),
 ('fc','laser.fc_hz','Laser cutoff frequency [Hz]','F1 high-frequency rolloff.','positive'),
 ('c4','laser.C4','Cavity C4 [rad² Hz³]','F3 1/f⁴ coefficient.','nonnegative'),
 ('c3','laser.C3','Cavity C3 [rad² Hz²]','F3 1/f³ coefficient.','nonnegative'),
 ('c2','laser.C2','Cavity C2 [rad² Hz]','F3 1/f² coefficient.','nonnegative'),
 ('bHz','laser.B_hz','Laser servo B [Hz]','F4 servo parameter; supply an apparatus value.','positive'),
 ('gamma','laser.gamma','Laser servo gamma [dimensionless]','F4 servo numerator parameter.','positive'),
 ('delta','laser.delta','Laser servo delta [dimensionless]','F4 servo denominator parameter.','positive'),
 ('fiberL','line.l','Fiber noise per length [rad² Hz/km]','Eq.6 noise amplitude per single-pass kilometer.','nonnegative'),
 ('fc1','line.fc1_hz','Fiber cutoff frequency [Hz]','Eq.6 high-frequency rolloff.','positive'),
 ('indexN','physics.n','Fiber refractive index [dimensionless]','Propagation index; blank uses the sourced reference value with a report notice.','positive'),
 ('kFactor','physics.K','Pass-correlation factor K [dimensionless]','Eq.5: 4 for correlated passes, 2 otherwise.','correlation'),
 ('lambdaS','physics.lambda_s_nm','Measurement wavelength [nm]','Dual-band wavelength in Eq.8.','positive'),
 ('lambdaQ','physics.lambda_q_nm','Quantum wavelength [nm]','Quantum-channel wavelength in Eq.8.','positive'),
 ('s0','physics.s0','Dual-band detection PSD amplitude [rad²/Hz]','Added once at the central node, not multiplied by K (Appendix G).','nonnegative'),
 ('fc2','physics.fc2_hz','Dual-band detection cutoff [Hz]','Detection noise rolloff in Eq.8.','positive'),
 ('omegaA','actuator.omega_a_rad_s','Actuator pole angular frequency [rad/s]','One-pole actuator model; an explicit engineering approximation.','positive'),
 ('actuatorBandwidth','actuator.bandwidth_hz','Actuator pole frequency [Hz]','Converted to omega_a; not the round-trip loop suppression bandwidth.','positive'),
 ('sigmaLimit','operation.sigma_limit_rad','Phase RMS threshold [rad]','Spectral operating-window threshold in Eq.4.','phase_positive'),
 ('tauMax','operation.tau_max_s','Maximum key window [s]','Upper cap on the solved operating window, not a fixed operating time.','positive'),
 ('operationTauPs','operation.tau_ps_s','Time spent on phase stabilization between key windows [s]','Downtime entering the duty cycle.','nonnegative'),
 ('lineCsvLength','line.spectrum.measurement_length_km','Measured fiber length [km]','Length associated with the uploaded single-arm PSD.','positive')]:field(id,path,label,help,rule)
for prefix in ('laser','line'):
    for suffix,label,help in [
     ('CsvFile','Spectrum CSV filename','Match an uploaded frequency/PSD CSV filename.'),
     ('CsvMode','Spectrum use','Fit the published model, or interpolate the measured PSD directly.'),
     ('CsvQuantity','PSD quantity','Identify phase or frequency noise; conversion requires this metadata.'),
     ('CsvPsdUnit','PSD units','Use rad^2/Hz for phase or Hz^2/Hz for frequency.'),
     ('CsvSidedness','PSD sidedness','One- or two-sided normalization is mandatory.'),
     ('CsvPass','Measurement pass','Single or round-trip; round-trip also needs an explicit PSD multiplier.'),
     ('RoundTripFactor','Round-trip PSD multiplier [dimensionless]','Measured round-trip PSD divided by this explicit factor gives the single-pass model.'),
     ('Extrapolation','Outside measured band','Reject by default, or explicitly allow endpoint power-law extrapolation.')]:
        field(prefix+suffix,None,label,help,'positive' if suffix=='RoundTripFactor' else None)
for id,label,help in [
 ('mode','Phase input mode','Measured phase bypasses spectra; spectral forecast integrates explicit source models.'),
 ('protocol','Key protocol','Choose an asymptotic protocol.'),('detectorMode','Detector representation','One scalar receiver or two separate scalar projections.'),
 ('reach','Optional reach scan','Off by default. A single observation or complete-loss point does not determine length dependence.'),
 ('laserModel','Laser model','F1 free laser or F2–F4 cavity-stabilized laser.'),('laserInput','Laser data source','Choose direct coefficients, linewidth or a spectrum with normalization metadata.'),
 ('schemeLasers','Laser arrangement','Eq.5 common laser or Eq.7 independent lasers.'),('schemeCompensation','Fiber compensation','Free, Eq.8 dual-band or Williams remote-end classical loop.'),
 ('actuatorInput','Actuator representation','Supply the pole, pole frequency, or a measured frequency response.'),
 ('actuatorCsvFile','Actuator response CSV filename','Upload frequency, magnitude and phase, not a noise PSD.'),
 ('actuatorMagnitudeUnit','Response magnitude unit','Linear ratio or decibels.'),('actuatorPhaseUnit','Response phase unit','Degrees for the one-pole response fit.'),
 ('extraToml','Additional TOML','Explicit ranges and numerical settings; avoid duplicating generated sections.'),
 ('files','CSV uploads','Files stay on the local server for this report session.'),
 ('toml','Configuration editor — core convention','Advanced input: SNS decoy intensities here are already two-user values; never double them again.'),
 ('labelA','Variant A name','A distinct short name for this calculation.'),('labelB','Variant B name','A distinct short name for this calculation.'),('labelC','Variant C name','Optional third calculation name.'),
 ('tomlA','Variant A TOML','Explicit core configuration; SNS intensities use the two-user convention.'),('tomlB','Variant B TOML','Explicit core configuration; SNS intensities use the two-user convention.'),('tomlC','Variant C TOML','Optional explicit core configuration; SNS intensities use the two-user convention.'),
 ('lengths','Working total lengths [km]','Comma-separated positive total lengths for the spectral comparison.'),
 ('exampleChoice','Published example','Load explicit source values; these are examples, not apparatus defaults.')]:field(id,None,label,help)

RULE_MESSAGES = dict(positive='Must be > 0.',nonnegative='Must be >= 0.',probability='Must be in (0, 1].',sns_probability='SNS requires both sending and not-sending events; epsilon must be strictly between 0 and 1. Published values are typically 0.05-0.3.',error='Must be in [0, 0.5).',phase='Must be in [0, pi] rad.',phase_positive='Must be in (0, pi] rad.',fec='Must be >= 1.',integer='Must be a nonnegative integer.',correlation='Must be 2 or 4.',finite='Must be finite.')

def valid_number(number, rule):
    if isinstance(number,bool) or not isinstance(number,(int,float)) or not math.isfinite(number):return False
    return {'positive':lambda:number>0,'nonnegative':lambda:number>=0,'probability':lambda:0<number<=1,'sns_probability':lambda:0<number<1,
            'error':lambda:0<=number<.5,'phase':lambda:0<=number<=math.pi,'phase_positive':lambda:0<number<=math.pi,
            'fec':lambda:number>=1,'integer':lambda:number>=0 and int(number)==number,
            'correlation':lambda:number in (2,4),'finite':lambda:True}[rule]()

def unusual_message(number, kind):
    if kind=='high_phase' and number>1:return 'Unusually large phase RMS (> 1 rad); calculation is allowed.'
    if kind=='low_efficiency' and number<.05:return 'Unusually low detector efficiency (< 0.05); calculation is allowed.'
    if kind=='high_fec' and number>2:return 'Unusually large f_EC (> 2); calculation is allowed.'
    if kind=='sns_typical' and not .05<=number<=.3:return 'Outside the typical sending range 0.05-0.3; calculation is allowed.'
    return None

class FieldValidationError(ValueError):
    def __init__(self, errors):
        self.field_errors=errors
        super().__init__('Correct the highlighted input fields before calculating.')

def get_path(raw, path):
    item=raw
    try:
        for key in path.split('.'):
            item=item[int(key)] if isinstance(item,list) else item[key]
        return item
    except (KeyError,IndexError,TypeError,ValueError):return None

def validate_web_input(raw):
    errors={};warnings={}
    def check(path,number,rule,unusual=None):
        if number is None:return
        if not valid_number(number,rule):errors[path]=RULE_MESSAGES[rule];return
        message=unusual_message(number,unusual)
        if message:warnings[path]=message
    for spec in FIELDS.values():
        if spec['path'] and spec['rule']:check(spec['path'],get_path(raw,spec['path']),spec['rule'],spec['unusual'])
    # Direct core-keyrate and physical representations are also accepted in Advanced.
    aliases={'detector_efficiency':('probability','low_efficiency'),'detector_dark_count_rate_hz':('nonnegative',None),
             'detector_error':('error',None),'attenuation_db_per_km':('positive',None)}
    for name,(rule,unusual) in aliases.items():check('keyrate.'+name,get_path(raw,'keyrate.'+name),rule,unusual)
    for spec in FIELDS.values():
        if spec['path'] and spec['rule'] and spec['path'].startswith(('laser.','line.')):
            name=spec['path'].split('.')[-1]
            if name in ('r3','r2','fc_hz','C4','C3','C2','B_hz','gamma','delta','l','fc1_hz'):
                check('physics.'+name,get_path(raw,'physics.'+name),spec['rule'])
    for name in ('arm_a_km','arm_b_km'):check('line.'+name,get_path(raw,'line.'+name),'positive')
    for node in ('laser','line'):
        check(node+'.spectrum.round_trip_psd_factor',get_path(raw,node+'.spectrum.round_trip_psd_factor'),'positive')
    p=raw.get('keyrate',{})
    decoys=[p.get(k) for k in ('decoy_big','decoy_medium','decoy_mini')]
    if all(isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x) for x in decoys):
        if not decoys[0]>decoys[1]>decoys[2]>=0:
            for key in ('decoy_big','decoy_medium','decoy_mini'):errors['keyrate.'+key]='Require big > medium > mini >= 0.'
    line=raw.get('line',{})
    for key,other in [('loss_a_db','loss_b_db'),('loss_b_db','loss_a_db')]:
        if key in line and other not in line:errors['line.'+other]='Supply both complete arm losses.'
    total,imbalance=line.get('length_km'),line.get('imbalance_km')
    if isinstance(total,(int,float)) and isinstance(imbalance,(int,float)) and abs(imbalance)>=total:
        errors['line.imbalance_km']='Magnitude must be less than total length (both arms positive).'
    a,b=p.get('nmin_cal'),p.get('nmax_cal')
    if isinstance(a,(int,float)) and isinstance(b,(int,float)) and a>=b:errors['keyrate.nmax_cal']='Must be greater than n_min.'
    if errors:raise FieldValidationError(errors)
    return warnings
