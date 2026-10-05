"""Explicit laboratory inputs; no equipment presets or inherited apparatus values.

Run: python -m tfqkd.lab_inputs INPUT.toml
Physical reference values belong to examples, not this resolver.
"""
from __future__ import annotations
import copy
from dataclasses import dataclass, field
import json
from pathlib import Path
import tomllib
import numpy as np
from .config import ROOT
from .measured_inputs import fit_psd, fit_actuator, MeasurementRejected

NONNEGATIVE_PHYSICS = frozenset(('r3','r2','C4','C3','C2','l','s0'))
MONOTONE_CANDIDATES = tuple('physics.'+key for key in sorted(NONNEGATIVE_PHYSICS)) + ('keyrate.detector_dark_count_rate_hz',)
CORNING_SOURCE = 'https://www.corning.com/content/dam/corning/media/worldwide/coc/documents/Fiber/product-information-sheets/PI-1424-AEN.pdf'
DEFAULT_N = 1.4682 # Corning SMF-28 Ultra PI-1424-AEN (July 2025), p.2, effective GROUP index at 1550 nm.
SPEED_OF_LIGHT_KM_S = 299792.458 # QKD.ipynb Cell 3; exact speed of light expressed in km/s.
DEFAULT_K = 4 # Eq.5 discussion, bertaina2024; correlated forward/backward phase noise, user-overridable.
PHYSICS_FIELDS = NONNEGATIVE_PHYSICS | {'fc_hz','B_hz','gamma','delta','fc1_hz','fc2_hz','lambda_s_nm','lambda_q_nm','c_km_s','n','K'}
KEYRATE_COMMON = {'clockrate_hz','attenuation_db_per_km','detector_efficiency','detector_dark_count_rate_hz','detector_error','f_error'}
KEYRATE_SNS = {'decoy_big','decoy_medium','decoy_mini','pz_sns','eps_sns_aopp'}
KEYRATE_CAL = {'pz_cal','nmin_cal','nmax_cal','u_cal'}
KEYRATE_FIELDS = KEYRATE_COMMON | KEYRATE_SNS | KEYRATE_CAL | {'stab_overhead_s','duty_bb84'}
OPERATION_FIELDS = {'sigma_limit_rad','tau_max_s','tau_ps_s'}

class MissingInputs(ValueError):
    def __init__(self, parameters):
        self.parameters=sorted(set(parameters))
        super().__init__('Missing explicit inputs: '+', '.join(self.parameters))

def read_toml(path):
    with Path(path).open('rb') as stream:
        return tomllib.load(stream)


def merge(base, extra):
    result = copy.deepcopy(base)
    for key, value in extra.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def lorentz_width_to_r2(width_hz):
    if not np.isfinite(width_hz) or width_hz < 0:
        raise ValueError('Lorentz white-noise FWHM must be finite and nonnegative [Hz]')
    # Eqs. 1,5, didomenico2010 (single-sided S_nu=f^2*S_phi, FWHM=pi*h0); Eq. F1, bertaina2024: h0=r2.
    return float(width_hz/np.pi)


@dataclass
class LabInput:
    config: dict
    settings: dict
    provenance: list
    ranges: dict
    warnings: list
    missing_amplitudes: tuple = ()
    direct_spectra: dict = field(default_factory=dict)

    def clone(self):
        return copy.deepcopy(self)

    def arms(self):
        line = self.settings['line']
        if 'arm_a_km' in line:
            return float(line['arm_a_km']), float(line['arm_b_km'])
        # Eq. 5, bertaina2024, definition delta_L=L_A-L_B; geometry identity for a physical total length.
        return ((line['length_km']+line['imbalance_km'])/2,
                (line['length_km']-line['imbalance_km'])/2)

    def set_parameter(self, path, value):
        section, name = path.split('.')
        if section in ('physics','keyrate','operation'):
            self.config[section][name] = value
        else:
            self.settings[section][name] = value
        if path == 'laser.lorentz_width_hz':
            self.config['physics']['r2'] = lorentz_width_to_r2(value)
        if path == 'operation.tau_ps_s':
            self.config['keyrate']['stab_overhead_s'] = value
        if path in self.missing_amplitudes:
            self.missing_amplitudes=tuple(p for p in self.missing_amplitudes if p!=path)
        if section == 'keyrate' and name in self.settings.get('_ideal_protocol_fields', []):
            self.settings['_ideal_protocol_fields'].remove(name)
        validate(self)


def validate(lab):
    if lab.settings.get('phase',{}).get('mode')=='measured':
        from .lab_phase_input import validate_phase
        return validate_phase(lab)
    c, s = lab.config, lab.settings
    a, b = lab.arms()
    if not np.all(np.isfinite([a,b])) or min(a,b) <= 0:
        raise ValueError('Both physical arm lengths must be finite and positive; |imbalance| < total length')
    for name, value in c['physics'].items():
        if not np.isfinite(value) or (value < 0 if name in NONNEGATIVE_PHYSICS else value <= 0):
            raise ValueError(f'Invalid physical input {name}={value}')
    if c['physics']['K'] not in (2,4):
        raise ValueError('K must be 2 or 4 (Eq. 5 discussion, bertaina2024)')
    for name, value in c['keyrate'].items():
        if not np.isfinite(value) or value < 0:
            raise ValueError(f'Invalid key-rate input {name}')
    p = c['keyrate']
    if not 0 < p['detector_efficiency'] <= 1 or not 0 <= p['detector_error'] < 0.5:
        raise ValueError('Require 0 < detector efficiency <= 1 and 0 <= detection error < 0.5')
    if not 0 < p['clockrate_hz'] or p['detector_dark_count_rate_hz'] >= p['clockrate_hz']:
        raise ValueError('Clock must be positive and dark count probability below one')
    if (s['protocol']['name']=='SNS-AOPP' and not p['decoy_big'] > p['decoy_medium'] > p['decoy_mini'] >= 0) or p['f_error'] < 1:
        raise ValueError('Require ordered decoys and error correction factor >= 1')
    if (s['protocol']['name']=='SNS-AOPP' and not 0 < p['eps_sns_aopp'] < 1) or not 0 <= p['pz_sns'] <= 1 or not 0 <= p['pz_cal'] <= 1:
        raise ValueError('Invalid protocol probabilities')
    if s['protocol']['name']=='CAL' and (not isinstance(p['nmin_cal'],int) or not isinstance(p['nmax_cal'],int) or not 0 <= p['nmin_cal'] < p['nmax_cal']):
        raise ValueError('CAL truncation must use ordered nonnegative integers')
    if s['protocol']['name']=='CAL' and p['u_cal']<=0:
        raise ValueError('CAL intensity u_cal must be positive')
    op = c['operation']
    if not all(np.isfinite(op[key]) for key in ('sigma_limit_rad','tau_max_s','tau_ps_s')) or op['sigma_limit_rad'] <= 0 or op['tau_max_s'] <= 0 or op['tau_ps_s'] < 0:
        raise ValueError('Invalid phase threshold or operation times')
    if 1/op['tau_max_s'] < c['grid']['f_min_hz']:
        raise ValueError('Time cap exceeds the supported lower frequency bound')
    if s['scheme']['lasers'] not in ('common','independent') or s['scheme']['compensation'] not in ('none','dual','classical'):
        raise ValueError('Unknown laser arrangement or compensation scheme')
    if s['protocol']['name'] not in ('SNS-AOPP','CAL'):
        raise ValueError('Protocol must be SNS-AOPP or CAL')
    k = s['loop']['safety_fraction']
    if not np.isfinite(k) or not 0 < k < 1:
        raise ValueError('Automatic-g safety_fraction must lie strictly between zero and one')
    if 'g_per_s' in s['loop'] and (not np.isfinite(s['loop']['g_per_s']) or s['loop']['g_per_s'] < 0):
        raise ValueError('Manual g must be finite and nonnegative; zero disables feedback')
    from .lab_observed import explicit_losses,receiver_inputs,channel_noise
    explicit_losses(s['line'])
    channels=receiver_inputs(s['detector'])
    if channels and any(channel_noise(ch)>=p['clockrate_hz'] for ch in channels):raise ValueError('Each noise-count probability must be below one')
    if s['scheme']['compensation'] == 'classical' and 'line' not in lab.direct_spectra:
        omega = s['actuator'].get('omega_a_rad_s')
        if omega is None or not np.isfinite(omega) or omega <= 0:
            raise ValueError('Classical compensation needs a positive explicit omega_a_rad_s or measured response')


def resolve(user=None, *, core_config=None):
    """Resolve only explicit physical inputs. core_config supplies numerical settings ONLY."""
    from .keyrates import KeyRateParameters
    raw=read_toml(user) if isinstance(user,(str,Path)) else copy.deepcopy(user or {})
    base_dir=Path(user).resolve().parent if isinstance(user,(str,Path)) else Path.cwd()
    if 'phase' in raw:
        from .lab_phase_input import resolve_phase
        return resolve_phase(raw,core_config=core_config)
    defaults=read_toml(ROOT/'configs/lab_defaults.toml')
    allowed=set(defaults)|{'line','detector','actuator','physics','keyrate','operation','ranges'}
    if set(raw)-allowed:raise ValueError(f'Unknown input sections: {sorted(set(raw)-allowed)}')
    s=merge(defaults,{k:v for k,v in raw.items() if k not in ('physics','keyrate','operation','ranges')})
    for section in ('line','detector','actuator'):s.setdefault(section,{})
    fields=dict(laser={'model','lorentz_width_hz','spectrum'}|PHYSICS_FIELDS.intersection({'r3','r2','fc_hz','C4','C3','C2','B_hz','gamma','delta'}),
                line={'length_km','imbalance_km','arm_a_km','arm_b_km','loss_a_db','loss_b_db','l','fc1_hz','attenuation_db_per_km','spectrum'},
                detector={'efficiency','dark_count_rate_hz','error','channels'},
                actuator={'omega_a_rad_s','bandwidth_hz','response'},scheme={'lasers','compensation'},
                protocol={'name'},loop={'safety_fraction','g_per_s'},operation=OPERATION_FIELDS,
                physics=PHYSICS_FIELDS,keyrate=KEYRATE_FIELDS)
    for section,keys in fields.items():
        unknown=set(raw.get(section,{}))-keys
        if unknown:raise ValueError(f'Unknown fields in {section}: {sorted(unknown)}; no equipment preset selectors are supported')
    c={k:copy.deepcopy(s[k]) for k in ('grid','performance','validation','classical','classical_actuator')}
    if core_config is not None:
        for section in c:c[section]=merge(c[section],core_config.get(section,{}))
        for section in c:c[section]=merge(c[section],raw.get(section,{}))
        if 'scenarios' in core_config:c['scenarios']=copy.deepcopy(core_config['scenarios'])
    c.update(physics=copy.deepcopy(raw.get('physics',{})),keyrate=copy.deepcopy(raw.get('keyrate',{})),operation=copy.deepcopy(raw.get('operation',{})))
    provenance=[];messages=[]
    from .lab_observed import receiver_inputs,explicit_losses,direct_psd,channel_noise
    channels=receiver_inputs(s['detector']);losses=explicit_losses(s['line']);direct={}
    def put(section,key,value,label):
        if key in c[section]:raise ValueError(f'Duplicate representation for {section}.{key}; choose exactly one input method')
        c[section][key]=value
        provenance.append(dict(parameter=section+'.'+key,value=value,input=label))
    for section in ('physics','keyrate','operation'):
        for key,value in c[section].items():provenance.append(dict(parameter=section+'.'+key,value=value,input='explicit input; source not inferred'))
    for key in fields['laser']-{'model','lorentz_width_hz','spectrum'}:
        if key in raw.get('laser',{}):put('physics',key,s['laser'][key],'direct laser coefficient')
    for key in ('l','fc1_hz'):
        if key in s['line']:put('physics',key,s['line'][key],'direct line coefficient')
    if 'attenuation_db_per_km' in s['line']:put('keyrate','attenuation_db_per_km',s['line']['attenuation_db_per_km'],'explicit line attenuation')
    for key,target in [('efficiency','detector_efficiency'),('dark_count_rate_hz','detector_dark_count_rate_hz'),('error','detector_error')]:
        if key in s['detector']:put('keyrate',target,s['detector'][key],'explicit detector characteristic')
    if channels:
        put('keyrate','detector_efficiency',channels[0]['efficiency'],'D0 explicit storage; both detector projections calculated separately')
        put('keyrate','detector_dark_count_rate_hz',channel_noise(channels[0]),'D0 explicit dark plus supplied background storage; not a receiver average')
        messages.append('SCALAR_DETECTOR_PROJECTIONS: two real channels retained; two scalar receiver projections, no averaging and no certified bounds for unequal detectors.')
    if losses:
        if 'attenuation_db_per_km' in c['keyrate']:raise ValueError('Complete arm losses conflict with per-km attenuation')
        c['keyrate']['attenuation_db_per_km']=0. # Unused storage: explicit TOTAL losses are passed to the unchanged protocol.
        provenance.append(dict(parameter='line.complete_arm_losses_db',value=losses,input='explicit losses including insertions; author equalization to lossier arm'))
    if 'lorentz_width_hz' in s['laser']:
        put('physics','r2',lorentz_width_to_r2(s['laser']['lorentz_width_hz']),'pure Lorentz component; Eqs.1,5 didomenico2010 / Eq.F1 bertaina2024')
        messages.append('Lorentz width supplies only r2 of the free laser plant; no r3 or fc is inferred.')
    for node in ('laser','line'):
        if 'spectrum' in s[node]:
            if s[node]['spectrum'].get('mode','fit')=='direct':
                data,metadata=direct_psd(s[node]['spectrum'],node,base_dir,c['grid'])
                direct[node]=data;provenance.append(dict(parameter=node+'.spectrum',**metadata))
                messages.append(f'DIRECT_PSD: {node} supplied without parametric fit; inverse requirements and coefficient sensitivity are not calculated for this node.')
                if metadata['extrapolation_used']:
                    messages.append(f'DIRECT_PSD_EXTRAPOLATION: {node} extrapolated outside {metadata["measured_band_hz"]} Hz by endpoint power laws {metadata["endpoint_log_slopes"]}; engineering assumption, not measured data.')
                if node=='line':messages.append('DIRECT_ARM_PSD: input is the fiber-only arm PSD for the SELECTED compensation state. Other schemes and loop settings cannot be predicted from it. Linear length scaling follows the Eq.6/8 model assumption, not a calibrated spatial map.')
                continue
            if s[node]['spectrum'].get('mode','fit')!='fit':raise MeasurementRejected('PSD mode must be fit or direct')
            if node=='laser' and s['laser']['model']!='free':raise MeasurementRejected('F1 fitting describes a FREE laser plant; do not fit a locked output as its free plant')
            parameters,metadata=fit_psd(s[node]['spectrum'],node,s['fit'],base_dir)
            for key,value in parameters.items():put('physics',key,value,'measured CSV fit')
            provenance.append(dict(parameter=node+'.spectrum',**metadata))
            messages.append('B6_MODEL_RESIDUAL: CSV fit is a model approximation; conditional fit SE is not a measured equipment uncertainty.')
    actuator=s['actuator']
    representations=set(actuator)&{'omega_a_rad_s','bandwidth_hz','response'}
    if len(representations)>1:raise ValueError('Choose actuator omega_a_rad_s OR bandwidth_hz OR response CSV')
    if 'bandwidth_hz' in actuator:
        if not np.isfinite(actuator['bandwidth_hz']) or actuator['bandwidth_hz']<=0:raise ValueError('Actuator bandwidth must be positive [Hz]')
        # User-adopted engineering pole; conversion omega_a=2*pi*f_a is a unit identity, not an apparatus identification.
        actuator['omega_a_rad_s']=float(2*np.pi*actuator['bandwidth_hz'])
        provenance.append(dict(parameter='actuator.omega_a_rad_s',value=actuator['omega_a_rad_s'],input='explicit nominal bandwidth; engineering one-pole proxy'))
    if 'response' in actuator:
        omega,metadata=fit_actuator(actuator['response'],s['fit'],base_dir)
        actuator['omega_a_rad_s']=omega;provenance.append(dict(parameter='actuator.response',**metadata))
    if s['laser']['model'] not in ('free','cavity'):raise ValueError('laser.model must be free or cavity')
    missing=[];line=s['line']
    if {'arm_a_km','arm_b_km'} & set(line):
        if {'length_km','imbalance_km'} & set(line):raise ValueError('Use both arms OR total length plus imbalance, without conflicting geometry')
        missing.extend('line.'+key for key in ('arm_a_km','arm_b_km') if key not in line)
    else:
        missing.extend('line.'+key for key in ('length_km','imbalance_km') if key not in line)
    required_physics={'fc_hz','fc1_hz'}
    amplitude_names={'r3','r2','l'}
    if s['laser']['model']=='cavity':
        required_physics|={'B_hz','gamma','delta'};amplitude_names|={'C4','C3','C2'}
    if s['scheme']['compensation']=='dual':
        required_physics|={'lambda_s_nm','lambda_q_nm','fc2_hz'};amplitude_names.add('s0')
    if 'laser' in direct:
        required_physics-={'fc_hz','B_hz','gamma','delta'};amplitude_names-={'r3','r2','C4','C3','C2'}
    if 'line' in direct:
        required_physics-={'fc1_hz','lambda_s_nm','lambda_q_nm'};amplitude_names.discard('l')
        if 'g_per_s' in s['loop']:raise ValueError('Direct arm PSD describes a fixed measured state; remove the manual g request because another controller cannot be inferred')
    missing.extend('physics.'+key for key in required_physics if key not in c['physics'])
    from .lab_protocol_requirements import ideal_protocol_inputs, validate_target
    ideal_protocol_inputs(c['keyrate'], s, provenance, messages)
    validate_target(s)
    required_keyrate=(KEYRATE_COMMON-{'detector_dark_count_rate_hz'})|(KEYRATE_SNS if s['protocol']['name']=='SNS-AOPP' else KEYRATE_CAL)
    missing.extend('keyrate.'+key for key in required_keyrate if key not in c['keyrate'])
    missing.extend('operation.'+key for key in OPERATION_FIELDS if key not in c['operation'])
    if s['scheme']['compensation']=='classical' and 'line' not in direct and 'omega_a_rad_s' not in actuator:missing.append('actuator.omega_a_rad_s (or bandwidth_hz/response)')
    if missing:raise MissingInputs(missing)
    unknown_amplitudes=[];ranges=copy.deepcopy(raw.get('ranges',{}))
    for section,names in [('physics',amplitude_names),('keyrate',{'detector_dark_count_rate_hz'})]:
        for key in sorted(names):
            path=section+'.'+key
            if key not in c[section]:
                if path in ranges:c[section][key]=ranges[path]['minimum']
                else:
                    c[section][key]=0.;unknown_amplitudes.append(path)
                    provenance.append(dict(parameter=path,value=0.,input='UNMEASURED: optimistic zero contribution, not a measured value or default'))
    if unknown_amplitudes:
        messages.append('UNMEASURED_AMPLITUDES: '+', '.join(unknown_amplitudes)+'; zero contributions give an upper key estimate, with separate conditional requirements.')
    if 'n' not in c['physics']:
        c['physics']['n']=DEFAULT_N
        provenance.append(dict(parameter='physics.n',value=DEFAULT_N,default_used=True,source=CORNING_SOURCE,page=2,
                               reference='Effective group index at 1550 nm; Corning SMF-28 Ultra PI-1424-AEN, July 2025'))
        messages.append('N_REFERENCE_DEFAULT: n=1.4682 (group index at 1550 nm, Corning p.2). Archived Table I sensitivity for n=1.44–1.48 reaches 1.818719% in variance relative to n=1.45; this is not a bound for arbitrary user geometry.')
    c['physics'].setdefault('c_km_s',SPEED_OF_LIGHT_KM_S);c['physics'].setdefault('K',DEFAULT_K)
    # Unused dataclass fields are neutral storage only; the selected protocol never reads these placeholders.
    # Appendix A/B, bertaina2024: required protocol fields above remain explicitly supplied.
    for key in KEYRATE_FIELDS-required_keyrate-{'detector_dark_count_rate_hz','stab_overhead_s'}:c['keyrate'].setdefault(key,0)
    overhead=c['operation']['tau_ps_s']
    if 'stab_overhead_s' in c['keyrate'] and c['keyrate']['stab_overhead_s']!=overhead:raise ValueError('Conflicting overhead; use operation.tau_ps_s')
    c['keyrate']['stab_overhead_s']=overhead
    if set(c['keyrate'])!=set(KeyRateParameters.__dataclass_fields__):raise ValueError('Unexpected protocol field mapping')
    req=s['requirements']
    if not np.isfinite(req['loss_fraction']) or not 0<req['loss_fraction']<1:raise ValueError('requirements.loss_fraction must be in (0,1)')
    if not np.isfinite(req['root_rtol']) or req['root_rtol']<=0 or not isinstance(req['max_bracket_doublings'],int) or req['max_bracket_doublings']<1:
        raise ValueError('Invalid numerical inverse-search settings')
    num=s['numerics']
    if any(not isinstance(num[key],int) or num[key]<minimum for key,minimum in [('gain_samples',3),('gain_refinements',1),('workers',1)]):
        raise ValueError('Invalid numerical gain sample/refinement/worker count')
    if any(not np.isfinite(num[key]) or num[key]<=0 for key in ('gain_xatol_fraction','gain_rate_rtol')):
        raise ValueError('Invalid numerical gain tolerances')
    approved=set(defaults['sensitivity']['directions'])
    directions=s['sensitivity']['directions']
    if set(directions)-approved:raise ValueError('Sensitivity direction is not approved for: '+', '.join(sorted(set(directions)-approved)))
    for path,factor in directions.items():
        if factor not in (0.5,2.0,'halve_shortfall') or (factor=='halve_shortfall' and path!='keyrate.detector_efficiency'):
            raise ValueError(f'Unsupported sensitivity transformation for {path}')
    lab=LabInput(c,s,provenance,ranges,messages,tuple(unknown_amplitudes),direct);validate(lab)
    if c['grid']['mode'] not in ('fast','reference'):raise ValueError('Grid mode must be fast or reference')
    if not 0<c['grid']['f_min_hz']<c['grid']['f_max_hz']:raise ValueError('Frequency bounds must be positive and increasing')
    allowed_ranges={section+'.'+key for section in ('physics','keyrate','operation') for key in c[section]}
    allowed_ranges|={'line.length_km','line.imbalance_km','line.arm_a_km','line.arm_b_km','actuator.omega_a_rad_s','laser.lorentz_width_hz'}
    for path,interval in ranges.items():
        if path not in allowed_ranges:raise ValueError(f'Unknown range parameter {path}')
        if set(interval)!={'minimum','maximum','sources'}:raise ValueError(f'Range {path} requires minimum/maximum/sources')
        lo,hi=interval['minimum'],interval['maximum']
        if not np.all(np.isfinite([lo,hi])) or hi<lo or not isinstance(interval['sources'],list) or not interval['sources'] or any(not isinstance(v,str) or not v.strip() for v in interval['sources']):raise ValueError(f'Invalid or unsourced interval {path}')
        if path in ('operation.LB_km','keyrate.stab_overhead_s'):raise ValueError('Use laboratory geometry and operation.tau_ps_s')
        if path.startswith('line.') and (('arm_' in path) != ('arm_a_km' in line)):raise ValueError('Range geometry conflicts with input representation')
        for value in (lo,hi):
            probe=lab.clone();probe.set_parameter(path,value)
    return lab


def reference_input(config, scenario=None, *, protocol='SNS-AOPP', overrides=None):
    """Validation fixture ONLY: make every physical reference value explicit."""
    raw=read_toml(ROOT/'examples/bertaina2024_table3.toml')
    raw['physics']=copy.deepcopy(config['physics']);raw['keyrate']=copy.deepcopy(config['keyrate'])
    raw['operation']={key:config['operation'][key] for key in OPERATION_FIELDS}
    if scenario is not None:
        raw['line']=dict(arm_a_km=config['operation']['LB_km']+scenario['delta_L_km'],arm_b_km=config['operation']['LB_km'])
        raw['laser']=dict(model='cavity' if scenario['cavity'] else 'free')
        raw['scheme']=dict(lasers='common' if scenario['common'] else 'independent',compensation='dual' if scenario['stabilized'] else 'none')
    raw['protocol']=dict(name=protocol)
    if overrides:
        for section in ('line','laser','detector','actuator'): 
            if section in overrides:raw[section]=copy.deepcopy(overrides[section])
        raw=merge(raw,overrides)
    return resolve(raw,core_config=config)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path)
    args=parser.parse_args()
    try:
        lab=resolve(args.config)
        print(json.dumps(dict(settings=lab.settings,physics=lab.config['physics'],keyrate=lab.config['keyrate'],
                              provenance=lab.provenance,missing_amplitudes=lab.missing_amplitudes,warnings=lab.warnings),indent=2))
    except (MissingInputs,MeasurementRejected) as error:
        print(json.dumps(dict(error=str(error),missing_inputs=getattr(error,'parameters',None),diagnostics=getattr(error,'diagnostics',None)),indent=2));raise SystemExit(2)
