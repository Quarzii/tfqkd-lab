"""Reproducible numerical validation. --quick: unit suite; --full: unchanged T1-T7 and Figure 3."""

# CLI import bootstrap; no calculation settings are changed.
import sys as _sys
from pathlib import Path as _Path
_PROJECT_ROOT = _Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_PROJECT_ROOT))
ROOT = _PROJECT_ROOT
import argparse
import copy
import json
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter
import numpy as np
import scipy
from tfqkd.config import load
from tfqkd.engine import independent_jobs, scan
from tfqkd.integration import frequency_grid, PhaseIntegral
from tfqkd.spectra import components
from tfqkd.keyrates import KeyRateParameters, sns_aopp_rate_per_pulse, cal_rate_per_pulse
from tfqkd.protocol import duty_cycle

ROOT=_PROJECT_ROOT
OUT=ROOT/'results/stage3_a'


def scenario_calibration(job):
    c, scenario = job
    operation=c['operation'];params=KeyRateParameters(**c['keyrate'])
    reference_f=frequency_grid(c,mode='reference')
    reference=PhaseIntegral(reference_f,components(reference_f,scenario,c)['total'])
    tau,_=reference.operating_time(operation['sigma_limit_rad'],operation['tau_max_s'])
    times=np.r_[c['validation']['comparison_times_s'],tau]
    # Eq. 4, bertaina2024: compare both fixed-time integrals and the solved threshold time.
    reference_v=reference.variance(times)
    sigma=np.sqrt(reference.variance(tau))
    # Eq. 5 geometry / notebook calc_sigma_tau_loss; same fixed LB as the Table I reproduction.
    loss=(2*operation['LB_km']+scenario['delta_L_km'])*params.attenuation_db_per_km
    # Appendix A A1 / Appendix B B1 and unnumbered duty formula, bertaina2024.
    ref_rates=np.array([model(loss,sigma,params)*duty_cycle(tau,operation['tau_ps_s'])*params.clockrate_hz
                       for model in (sns_aopp_rate_per_pulse,cal_rate_per_pulse)])
    rows=[]
    for count in c['performance']['grid_candidates']:
        f=frequency_grid(c,points=count)
        integral=PhaseIntegral(f,components(f,scenario,c)['total'])
        actual_tau,status=integral.operating_time(operation['sigma_limit_rad'],operation['tau_max_s'])
        sigma=np.sqrt(integral.variance(actual_tau))
        rates=np.array([model(loss,sigma,params)*duty_cycle(actual_tau,operation['tau_ps_s'])*params.clockrate_hz
                        for model in (sns_aopp_rate_per_pulse,cal_rate_per_pulse)])
        # Numerical relative errors only; not physical formulas or fitted coefficients.
        errors=np.abs(integral.variance(times)/reference_v-1)
        rows.append(dict(points=count,variance_max_relative=float(np.max(errors)),tau_relative=float(abs(actual_tau/tau-1)),
                         key_max_relative=float(np.max(abs(rates/ref_rates-1))),tau_q_s=float(actual_tau),
                         variance_rad2=float(integral.variance(actual_tau)),key_bps=rates.tolist(),status=status))
    return dict(scenario=scenario['name'],reference_tau_s=float(tau),reference_variances_rad2=reference_v.tolist(),
                times_s=times.tolist(),reference_keys_bps=ref_rates.tolist(),rows=rows)


def calibration():
    c=load();OUT.mkdir(parents=True,exist_ok=True)
    values=independent_jobs(scenario_calibration,[(c,s) for s in c['scenarios']],c['performance']['workers'])
    summary=[]
    for index,count in enumerate(c['performance']['grid_candidates']):
        variance=max(v['rows'][index]['variance_max_relative'] for v in values)
        tau=max(v['rows'][index]['tau_relative'] for v in values)
        key=max(v['rows'][index]['key_max_relative'] for v in values)
        passed=(variance<=c['performance']['regression_variance_rtol'] and tau<=c['performance']['regression_time_rtol'] and key<=c['performance']['regression_key_rtol'])
        summary.append(dict(points=count,variance_max_relative=variance,tau_max_relative=tau,key_max_relative=key,passed=passed))
    fast=next(r['points'] for r in summary if r['passed'])
    evidence=dict(config=c,scope='Table I: fixed LB=100 km, five configured Figure 8 times plus reference threshold time; smallest passing nested 2^k+1 grid, not a universal minimum',
                  summary=summary,scenarios=values,selected_fast_points=fast,passed=fast==c['grid']['fast_points'])
    (OUT/'grid_regression.json').write_text(json.dumps(evidence,indent=2))
    print(json.dumps(summary,indent=2),flush=True)
    assert evidence['passed']
    # Scheduling identity on independent omega_a jobs; no parallelism inside a vectorized block.
    omegas=[c['classical_actuator']['omega_a_min_rad_s'],c['classical_actuator']['omega_a_max_rad_s']]
    a=scan(c,lengths=[c['classical']['length_km']],fractions=[c['classical_actuator']['gain_fractions'][0]],omegas=omegas,workers=1)
    b=scan(c,lengths=[c['classical']['length_km']],fractions=[c['classical_actuator']['gain_fractions'][0]],omegas=omegas,workers=c['performance']['workers'])
    for key in a:np.testing.assert_array_equal(a[key],b[key])
    (OUT/'process_identity.json').write_text(json.dumps(dict(passed=True,omega_jobs=len(omegas),workers=c['performance']['workers'],fields=list(a)),indent=2))


def full():
    out=OUT/'full';out.mkdir(parents=True,exist_ok=True)
    # Only execution mode changes. All physical parameters and former validation settings are retained.
    cfg=out/'config_reference.toml'
    cfg.write_text((ROOT/'configs/config.toml').read_text().replace('mode = "fast"','mode = "reference"'))
    commands=[('T1_T4',[sys.executable,'scripts/run/run_stage1.py','--config',str(cfg),'--output-dir',str(out)]),
              ('Figure3',[sys.executable,'scripts/run/run_keyrates.py','--config',str(cfg),'--outdir',str(out)]),
              ('author_keyrates',[sys.executable,'scripts/validate/validate_keyrates.py']),
              ('T6_T7',[sys.executable,'scripts/run/run_stage2_checks.py','--config',str(cfg),'--output-dir',str(out)]),
              ('T5',[sys.executable,'scripts/run/run_classical.py','--config',str(cfg),'--output-dir',str(out)])]
    checks=[]
    for name,command in commands:
        print('Full reference check:',name,flush=True)
        start=perf_counter()
        with (out/(name+'.log')).open('w') as stream:
            done=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
        check=dict(name=name,command=command,exit_code=done.returncode,seconds=perf_counter()-start)
        checks.append(check)
        (OUT/'full_suite.json').write_text(json.dumps(dict(grid_mode='reference',reference_points=load()['grid']['reference_points'],checks=checks),indent=2))
        print(check,flush=True)
        if done.returncode:raise RuntimeError(f'{name} failed; inspect {out/(name+".log")}')
    (out/'keyrate_reference.json').write_bytes((ROOT/'results/keyrate_reference.json').read_bytes())
    # Verify all former acceptance assertions, retaining T4 as the explicitly accepted validation limit.
    e=json.loads((out/'evidence.json').read_text());s=json.loads((out/'stage2_checks.json').read_text());t5=json.loads((out/'actuator_evidence.json').read_text())
    assert e['T1']['pass'] and e['T2']['revised']['pass'] and e['unit_tests_exit_code']==0
    assert all(r['pass'] for r in s['T6_white']) and s['T7_equal_arm_pass']
    assert max(t5[k] for k in ('max_accepted_time_change','max_accepted_cap_change','max_accepted_variance_change')) <= load()['classical']['quadrature_rtol']
    assert max(r['max_relative_error'] for r in t5['spatial_quadrature_checks']) <= load()['validation']['reference_rtol']
    checks.append(dict(name='former_acceptance_assertions',exit_code=0,T4='accepted discrepancy; unchanged validation limit'))
    (OUT/'full_suite.json').write_text(json.dumps(dict(grid_mode='reference',reference_points=load()['grid']['reference_points'],checks=checks,passed=True,
               versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)),indent=2))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--quick',action='store_true');parser.add_argument('--calibrate',action='store_true');parser.add_argument('--full',action='store_true');args=parser.parse_args()
    if args.quick:
        return subprocess.call([sys.executable,'-m','unittest','discover','-s','tests'],cwd=ROOT)
    if args.calibrate:calibration()
    if args.full:full()
    if not (args.quick or args.calibrate or args.full):parser.error('Choose --quick, --calibrate or --full')
    return 0

if __name__=='__main__':raise SystemExit(main())
