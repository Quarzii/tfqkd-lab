"""Compare local block 4 against unmodified notebook functions, identical inputs."""

# CLI import bootstrap; no calculation settings are changed.
import sys as _sys
from pathlib import Path as _Path
_PROJECT_ROOT = _Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_PROJECT_ROOT))
ROOT = _PROJECT_ROOT
import ast
import json
import math
from pathlib import Path
import numpy as np
from tfqkd.config import load
from tfqkd.keyrates import KeyRateParameters, cal_rate_per_pulse, sns_aopp_rate_per_pulse, decoy_bb84_rate_per_pulse

c=load(); p=KeyRateParameters(**c['keyrate'])
nb=json.loads((ROOT/'sources/data/QKD.ipynb').read_text()); namespace={'np':np,'math':math}
for index in (13,15,17,19,21):
    tree=ast.parse(''.join(nb['cells'][index]['source']))
    tree.body=[n for n in tree.body if isinstance(n,ast.FunctionDef)]
    exec(compile(tree,f'QKD.ipynb:cell{index}','exec'),namespace)
# QKD.ipynb cell 23 plot_panel_scenarios: parameter mapping.
shared=dict(pdc=p.detector_dark_count_rate_hz/p.clockrate_hz,DetectorEfficiency=p.detector_efficiency,fError=p.f_error)
cal=dict(shared,n_min=p.nmin_cal,n_max=p.nmax_cal,pZ=p.pz_cal,theta=2*np.arcsin(np.sqrt(p.detector_error)),alfa2=p.u_cal)
sns=dict(shared,eps=p.eps_sns_aopp,s=p.decoy_big/2,n=p.decoy_mini/2,pZ=p.pz_sns,decoy2=p.decoy_big,decoy1=p.decoy_medium,decoy0=p.decoy_mini,edet=p.detector_error)
bb=dict(shared,u=p.decoy_big,v=p.decoy_medium,w=p.decoy_mini,edet=p.detector_error,duty=p.duty_bb84)
rows=[]
for loss in np.linspace(c['keyrate_validation']['loss_min_db'],c['keyrate_validation']['loss_max_db'],c['keyrate_validation']['loss_points']):
    for sigma in (0.,c['operation']['sigma_limit_rad']):
        for name,own,ref,kwargs in [('CAL',cal_rate_per_pulse,'RdB_CAL_TF_QKD_err',cal),('SNS',sns_aopp_rate_per_pulse,'RdB_SNS_AOPP_TF_QKD_err',sns),('BB84',decoy_bb84_rate_per_pulse,'RdB_Decoy_err',bb)]:
            a=float(own(loss,sigma,p)); b=float(namespace[ref](loss,sigma,**kwargs))
            # Numerical comparison to QKD.ipynb functions, not a new physical equation.
            rows.append(dict(protocol=name,loss_db=loss,sigma=sigma,local=a,author=b,absolute_difference=abs(a-b),relative_difference=abs(a/b-1) if b else None))
summary={name:dict(max_abs=max(r['absolute_difference'] for r in rows if r['protocol']==name),max_relative=max(r['relative_difference'] for r in rows if r['protocol']==name and r['relative_difference'] is not None)) for name in ('CAL','SNS','BB84')}
(ROOT/'results/keyrate_reference.json').write_text(json.dumps(dict(summary=summary,rows=rows),indent=2))
print(json.dumps(summary,indent=2))
assert all(np.isclose(r['local'],r['author'],rtol=c['keyrate_validation']['relative_tolerance'],atol=np.finfo(float).eps) for r in rows)
