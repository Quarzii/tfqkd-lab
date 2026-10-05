"""User measurements, explicit normalization and adopted fit acceptance gates.

PSD fits reuse the B6 algorithm unchanged. Actuator fitting is the user's
engineering assumption, not a transfer function identified by Williams.
Run: python -m tfqkd.measured_inputs SPECIFICATION.json --node laser|line|actuator
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from .published_fit import fit_points

class MeasurementRejected(ValueError):
    def __init__(self,message,diagnostics=None):
        super().__init__(message);self.diagnostics=diagnostics or {}

def _table(spec,base_dir,columns):
    path=Path(spec['file']);path=path if path.is_absolute() else Path(base_dir)/path
    data=np.genfromtxt(path,delimiter=',',names=True,dtype=float,encoding='utf8')
    if data.dtype.names is None or not set(columns)<=set(data.dtype.names):
        raise MeasurementRejected(f'{path}: required CSV columns are {columns}')
    arrays=[np.atleast_1d(data[key]) for key in columns]
    if any(np.any(~np.isfinite(a)) for a in arrays):raise MeasurementRejected('CSV contains nonfinite values')
    order=np.argsort(arrays[0]);arrays=[a[order] for a in arrays]
    if len(arrays[0])<2 or np.any(arrays[0]<=0) or np.any(np.diff(arrays[0])<=0):
        raise MeasurementRejected('CSV needs distinct positive Fourier frequencies; omit DC and negative-frequency bins')
    return path,arrays

def normalize_psd(spec,base_dir='.'):
    required={'file','quantity','frequency_unit','psd_unit','sidedness','pass'}
    if required-set(spec):raise MeasurementRejected(f'PSD metadata missing: {sorted(required-set(spec))}')
    if spec['frequency_unit']!='Hz':raise MeasurementRejected('PSD frequency_unit must be Hz')
    quantity=spec['quantity']
    units={'phase':{'rad^2/Hz','rad2/Hz'},'frequency':{'Hz^2/Hz','Hz2/Hz'}}
    if quantity not in units or spec['psd_unit'] not in units[quantity]:
        raise MeasurementRejected('Specify phase PSD in rad^2/Hz or frequency PSD in Hz^2/Hz (linear density)')
    if spec['sidedness'] not in ('one-sided','two-sided') or spec['pass'] not in ('single','round-trip'):
        raise MeasurementRejected('Explicit sidedness one-sided/two-sided and pass single/round-trip are required')
    path,(f,psd)=_table(spec,base_dir,('frequency','psd'))
    if np.any(psd<=0):raise MeasurementRejected('Log PSD fit requires strictly positive linear density')
    factor=1.
    if spec['sidedness']=='two-sided':
        # Eq. B.1 and accompanying definition, li2013; Eqs.1,5, didomenico2010: one-sided real-noise density doubles.
        factor*=2
    if spec['pass']=='round-trip':
        if 'round_trip_psd_factor' not in spec:
            raise MeasurementRejected('Round-trip metadata also needs calibrated round_trip_psd_factor: correlated and uncorrelated paths are different (Eq.5, bertaina2024)')
        scale=spec['round_trip_psd_factor']
        if not np.isfinite(scale) or scale<=0:raise MeasurementRejected('round_trip_psd_factor must be positive')
        # Eq.5 discussion / Appendix G, bertaina2024, and Sec.4, clivati2018: divide by explicit calibrated PSD factor.
        # This factor is supplied by the user; four is NOT inferred for arbitrary geometry/frequency bands.
        factor/=scale
    elif 'round_trip_psd_factor' in spec:
        raise MeasurementRejected('round_trip_psd_factor conflicts with single-pass metadata')
    phase=psd*factor
    if quantity=='frequency':
        # Eq.1, didomenico2010; Eq.F1, bertaina2024; Snigirev2023 Methods: S_nu=f^2*S_phi.
        phase=phase/f**2
    return f,phase,dict(file=str(path.resolve()),input_metadata=spec,linear_normalization_factor=factor,
                        output_quantity='phase',output_unit='rad^2/Hz',output_sidedness='one-sided',output_pass='single')

def fit_psd(spec,node,settings,base_dir='.'):
    if node not in ('laser','line'):raise MeasurementRejected('PSD input is supported only for laser and line')
    f,phase,metadata=normalize_psd(spec,base_dir)
    model=dict(model='laser' if node=='laser' else 'fiber')
    if node=='line':
        if 'measurement_length_km' not in spec or not np.isfinite(spec['measurement_length_km']) or spec['measurement_length_km']<=0:
            raise MeasurementRejected('Fiber PSD fit needs measurement_length_km of the measured SINGLE-pass path')
        model['length_km']=spec['measurement_length_km']
    parameter_count=3 if node=='laser' else 2
    if len(f)<=parameter_count:raise MeasurementRejected('Too few PSD points to identify coefficients and residual scatter')
    fit,_=fit_points(f,phase,model,settings)
    metadata['fit']=fit
    threshold=spec.get('maximum_rms_log10_residual',settings['maximum_rms_log10_residual'])
    if not np.isfinite(threshold) or threshold<=0:raise MeasurementRejected('Fit RMS threshold must be positive [dex]')
    metadata['acceptance_threshold_dex']=threshold
    reasons=[]
    if fit['jacobian_rank']<parameter_count or any(v is None for v in fit['conditional_fit_se'].values()):
        reasons.append('parameters are not identifiable at the finite-difference Jacobian resolution')
    if fit['rms_log10_residual']>threshold:
        reasons.append(f'RMS residual {fit["rms_log10_residual"]:.9g} dex exceeds {threshold:.9g} dex')
    if reasons:raise MeasurementRejected('; '.join(reasons),metadata)
    return fit['parameters'],metadata

def fit_actuator(spec,settings,base_dir='.'):
    for key,value in [('frequency_unit','Hz'),('phase_unit','deg')]:
        if spec.get(key)!=value:raise MeasurementRejected(f'Actuator response requires {key}={value}')
    if spec.get('magnitude_unit') not in ('linear','dB'):
        raise MeasurementRejected('Actuator response needs magnitude_unit linear or dB; this is NOT a PSD')
    path,(f,magnitude,phase_deg)=_table(spec,base_dir,('frequency','magnitude','phase'))
    if spec['magnitude_unit']=='dB':
        # User-authorized engineering frequency-response input: amplitude dB=20*log10(|A|), not PSD dB.
        magnitude=10**(magnitude/20)
    if np.any(magnitude<=0) or len(f)<=1:raise MeasurementRejected('Require positive response magnitudes and multiple points')
    # User-authorized one-pole actuator model; degree-to-radian coordinate conversion is a unit identity.
    measured=magnitude*np.exp(1j*np.deg2rad(phase_deg))
    # User-adopted A(s)=1/(1+s/omega_a): an exactly flat unit response requires omega_a -> infinity,
    # so no finite pole is identified; machine epsilon only allows coordinate roundoff.
    if np.all(abs(measured-1)<=np.finfo(float).eps):
        raise MeasurementRejected('Flat unit response does not identify a finite actuator pole; extend the measured frequency band')
    def response(q):
        # Eq.A6, williams2008 + user-adopted actuator pole A(s)=1/(1+s/omega_a), NOT a sourced apparatus fit.
        return 1/(1+2j*np.pi*f/np.exp(q[0]))
    def residual(q):
        quotient=response(q)/measured
        return np.r_[np.log(abs(quotient)),np.angle(quotient)]
    margin=settings['log_parameter_search_margin']
    start=np.log(2*np.pi*np.median(f))
    opt=least_squares(residual,[np.clip(start,-margin,margin)],bounds=(-margin,margin),max_nfev=settings['max_nfev'])
    if not opt.success:raise MeasurementRejected('One-pole response fit did not converge')
    quotient=response(opt.x)/measured
    magnitude_rms=float(np.sqrt(np.mean(np.log10(abs(quotient))**2)))
    phase_rms=float(np.sqrt(np.mean(np.rad2deg(np.angle(quotient))**2)))
    # B6 adopted conditional local SE, here for complex log response with equal dimensionless residual weights.
    information=float((opt.jac.T@opt.jac)[0,0]);omega=float(np.exp(opt.x[0]))
    identifiable=information>np.finfo(float).eps*len(f)
    se=omega*np.sqrt(2*opt.cost/(len(opt.fun)-1)/information) if identifiable else None
    result=dict(file=str(path.resolve()),omega_a_rad_s=omega,conditional_fit_se_rad_s=se,
                rms_log_magnitude_dex=magnitude_rms,rms_phase_deg=phase_rms,
                engineering_assumption='unit DC gain one pole; joint log-magnitude and wrapped phase in radians; equal dimensionless residual weights',
                uncertainty_label='conditional local fit SE, not apparatus calibration uncertainty')
    magnitude_limit=spec.get('maximum_rms_log_magnitude_dex',settings['maximum_rms_log10_residual'])
    # User-approved engineering fit gate, not a source-provided phase calibration tolerance.
    phase_limit=spec.get('maximum_rms_phase_deg',settings['maximum_rms_phase_deg'])
    result.update(acceptance_phase_limit_deg=phase_limit,
                  acceptance_phase_default_used='maximum_rms_phase_deg' not in spec,
                  phase_gate_origin='user-approved engineering criterion; not from a publication')
    if not np.all(np.isfinite([magnitude_limit,phase_limit])) or min(magnitude_limit,phase_limit)<=0:
        raise MeasurementRejected('Response residual thresholds must be finite and positive')
    if not identifiable or magnitude_rms>magnitude_limit or phase_rms>phase_limit:
        raise MeasurementRejected(f'Response is not described by one pole: magnitude RMS={magnitude_rms:.9g} dex; phase RMS={phase_rms:.9g} deg; limits={magnitude_limit}, {phase_limit}',result)
    return omega,result

if __name__=='__main__':
    import argparse
    from .lab_inputs import read_toml
    from .config import ROOT
    parser=argparse.ArgumentParser();parser.add_argument('specification',type=Path);parser.add_argument('--node',choices=('laser','line','actuator'),required=True)
    args=parser.parse_args();spec=json.loads(args.specification.read_text());settings=read_toml(ROOT/'configs/lab_defaults.toml')['fit']
    try:
        output=fit_actuator(spec,settings,args.specification.parent) if args.node=='actuator' else fit_psd(spec,args.node,settings,args.specification.parent)
        print(json.dumps(dict(accepted=True,result=output),indent=2))
    except MeasurementRejected as error:
        print(json.dumps(dict(accepted=False,error=str(error),diagnostics=error.diagnostics),indent=2));raise SystemExit(2)
