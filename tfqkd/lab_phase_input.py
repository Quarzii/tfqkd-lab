"""Resolve a residual phase observation without constructing spectral parameters."""
import copy
import numpy as np
from .config import ROOT
from .lab_inputs import (LabInput, read_toml, merge, MissingInputs, KEYRATE_COMMON,
                         KEYRATE_SNS, KEYRATE_CAL, KEYRATE_FIELDS)
from .lab_observed import explicit_losses, receiver_inputs, channel_noise


def resolve_phase(raw, *, core_config=None):
    defaults = read_toml(ROOT/'configs/lab_defaults.toml')
    unknown = set(raw)-set(defaults)-{'phase','line','detector','actuator','physics','keyrate','operation','ranges'}
    if unknown:raise ValueError(f'Unknown input sections: {sorted(unknown)}')
    s = merge(defaults, {k:v for k,v in raw.items() if k not in ('physics','keyrate','operation','ranges')})
    for section in ('line','detector','actuator'):s.setdefault(section,{})
    phase = s['phase']
    allowed = {'mode','sigma_phi_rad','tau_s','tau_ps_s','description','source'}
    if set(phase)-allowed:raise ValueError(f'Unknown measured phase fields: {sorted(set(phase)-allowed)}')
    if phase.get('mode') != 'measured':raise ValueError('phase.mode must be measured')
    if any('spectrum' in raw.get(node,{}) for node in ('laser','line')):
        raise ValueError('Choose measured phase OR spectral inputs, not both')
    if raw.get('ranges'):raise ValueError('Measured phase mode accepts explicit point observations; specify separate observation files for distinct phase windows')
    line=s['line'];detector=s['detector'];protocol=s['protocol']['name']
    if protocol not in ('SNS-AOPP','CAL'):raise ValueError('Protocol must be SNS-AOPP or CAL')
    if set(line)-{'length_km','imbalance_km','arm_a_km','arm_b_km','loss_a_db','loss_b_db','attenuation_db_per_km'}:
        raise ValueError('Measured phase line input supports geometry and losses only')
    if set(detector)-{'efficiency','dark_count_rate_hz','error','channels'}:
        raise ValueError('Measured phase detector input supports explicit characteristics only')
    if set(raw.get('keyrate',{}))-KEYRATE_FIELDS:raise ValueError('Unknown key-rate fields')
    if set(raw.get('operation',{}))-{'tau_ps_s'}:raise ValueError('Measured phase mode uses phase.tau_s, no threshold or time cap')
    if 'tau_ps_s' in raw.get('operation',{}):
        if 'tau_ps_s' in phase:raise ValueError('Supply tau_PS once, in phase or operation')
        phase['tau_ps_s']=raw['operation']['tau_ps_s']
    missing=['phase.'+k for k in ('sigma_phi_rad','tau_s','tau_ps_s') if k not in phase]
    losses=explicit_losses(line);channels=receiver_inputs(detector)
    p=copy.deepcopy(raw.get('keyrate',{}));provenance=[];warnings=[];ideal=[]
    for key,value in p.items():provenance.append(dict(parameter='keyrate.'+key,value=value,input='explicit input'))
    def put(key,value,origin):
        if key in p:raise ValueError(f'Duplicate representation for keyrate.{key}')
        p[key]=value;provenance.append(dict(parameter='keyrate.'+key,value=value,input=origin))
    for key,target in [('efficiency','detector_efficiency'),('dark_count_rate_hz','detector_dark_count_rate_hz'),('error','detector_error')]:
        if key in detector:put(target,detector[key],'explicit detector input')
    if channels:
        put('detector_efficiency',channels[0]['efficiency'],'D0 storage; both scalar projections evaluated')
        put('detector_dark_count_rate_hz',channel_noise(channels[0]),'D0 dark plus explicit background storage; no detector average')
        warnings.append('SCALAR_DETECTOR_PROJECTIONS: two real channels retained; two hypothetical scalar receivers, not certified bounds for unequal detectors.')
    if losses:
        if 'attenuation_db_per_km' in p:raise ValueError('Use total losses OR per-km attenuation')
        p['attenuation_db_per_km']=0. # Unused dataclass storage, not an assumed measured attenuation.
        provenance.append(dict(parameter='line.complete_arm_losses_db',value=losses,input='explicit total arm losses with insertions'))
    elif 'attenuation_db_per_km' in line:put('attenuation_db_per_km',line['attenuation_db_per_km'],'explicit per-km loss')
    geometry=bool({'arm_a_km','arm_b_km','length_km','imbalance_km'}&set(line))
    if losses is None or geometry:
        if {'arm_a_km','arm_b_km'}&set(line):
            if {'length_km','imbalance_km'}&set(line):raise ValueError('Use two arms OR total length and imbalance')
            missing+=['line.'+k for k in ('arm_a_km','arm_b_km') if k not in line]
        else:missing+=['line.'+k for k in ('length_km','imbalance_km') if k not in line]
    from .lab_protocol_requirements import ideal_protocol_inputs, validate_target
    ideal_protocol_inputs(p, s, provenance, warnings)
    ideal = s['_ideal_protocol_fields']
    validate_target(s)
    required=(KEYRATE_COMMON-{'detector_dark_count_rate_hz'})|(KEYRATE_SNS if protocol=='SNS-AOPP' else KEYRATE_CAL)
    missing+=['keyrate.'+k for k in required if k not in p]
    if missing:raise MissingInputs(missing)
    unknown_amplitudes=[]
    if 'detector_dark_count_rate_hz' not in p:
        p['detector_dark_count_rate_hz']=0.;unknown_amplitudes.append('keyrate.detector_dark_count_rate_hz')
        warnings.append('UNMEASURED_AMPLITUDES: dark count contribution zero; conditional upper model estimate, not a measured zero.')
    for key in KEYRATE_FIELDS-required-{'detector_dark_count_rate_hz','stab_overhead_s'}:p.setdefault(key,0)
    if 'stab_overhead_s' in p and p['stab_overhead_s']!=phase['tau_ps_s']:raise ValueError('Conflicting tau_PS')
    p['stab_overhead_s']=phase['tau_ps_s']
    c={k:copy.deepcopy(s[k]) for k in ('grid','performance','validation','classical','classical_actuator')}
    if core_config:
        for k in c:c[k]=merge(c[k],core_config.get(k,{}))
    c.update(physics={},keyrate=p,operation={'tau_ps_s':phase['tau_ps_s']})
    if raw.get('physics') or raw.get('laser') or raw.get('actuator'):
        warnings.append('MEASURED_PHASE: supplied spectral/actuator parameters are not used; no spectral prediction or loop tuning is performed.')
    s['_ideal_protocol_fields']=ideal
    for key in ('sigma_phi_rad','tau_s','tau_ps_s'):
        provenance.append(dict(parameter='phase.'+key,value=phase[key],input='user-supplied residual phase / operating window',source=phase.get('source')))
    lab=LabInput(c,s,provenance,{},warnings,tuple(unknown_amplitudes));validate_phase(lab)
    return lab


def validate_phase(lab):
    phase=lab.settings['phase'];p=lab.config['keyrate'];line=lab.settings['line']
    if not np.all(np.isfinite([phase['sigma_phi_rad'],phase['tau_s'],phase['tau_ps_s']])) or phase['sigma_phi_rad']<0 or phase['tau_s']<=0 or phase['tau_ps_s']<0:
        raise ValueError('Require finite sigma_phi >= 0 [rad], tau > 0 [s], tau_PS >= 0 [s]')
    if ('arm_a_km' in line or 'length_km' in line):
        arms=lab.arms()
        if not np.all(np.isfinite(arms)) or min(arms)<=0:raise ValueError('Both arm lengths must be finite and positive')
    explicit_losses(line);channels=receiver_inputs(lab.settings['detector'])
    if not all(np.isfinite(v) and v>=0 for v in p.values()):raise ValueError('Invalid protocol/detector value')
    if not 0<p['detector_efficiency']<=1 or not 0<=p['detector_error']<.5 or p['clockrate_hz']<=0 or p['f_error']<1:
        raise ValueError('Invalid efficiency, clock, intrinsic error or fEC')
    darks=[channel_noise(ch) for ch in channels] if channels else [p['detector_dark_count_rate_hz']]
    if max(darks)>=p['clockrate_hz']:raise ValueError('Dark-count probability must be below one')
    if lab.settings['protocol']['name']=='SNS-AOPP':
        if not p['decoy_big']>p['decoy_medium']>p['decoy_mini']>=0 or not 0<p['eps_sns_aopp']<1 or not 0<=p['pz_sns']<=1:
            raise ValueError('Invalid SNS intensities/probabilities; decoys are total two-user intensities')
    elif not 0<=p['pz_cal']<=1 or p['u_cal']<=0 or not all(isinstance(p[k],int) for k in ('nmin_cal','nmax_cal')) or not 0<=p['nmin_cal']<p['nmax_cal']:
        raise ValueError('Invalid CAL intensity/probability/truncation')
