"""Independent spectral audit. No tfqkd import; source parameters are explicit.

Log-frequency adaptive quadrature and bracketed roots are deliberately different
from the core's trapezoid grid/interpolation. Protocol functions are extracted
unaltered from QKD.ipynb, not reimplemented here.
"""

# CLI import bootstrap; no calculation settings are changed.
import sys as _sys
from pathlib import Path as _Path
_PROJECT_ROOT = _Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_PROJECT_ROOT))
ROOT = _PROJECT_ROOT
import argparse
import ast
import json
import math
from pathlib import Path
import tomllib
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# Numerical audit controls, not physical coefficients; comparison tolerance is user-prescribed.
QUADRATURE_RTOL = 1e-11
ROOT_XTOL_HZ = 1e-10
ROOT_RTOL = 1e-13
QUADRATURE_SUBDIVISIONS = 1000


def spectra(f, p, a, b, *, common=False, stabilized=True, cavity=True):
    # Eq. F1, bertaina2024.
    free = p['r3']/f**3 + p['r2']/f**2*(p['fc_hz']/(f+p['fc_hz']))**2
    # Eq. F3, bertaina2024.
    cav = p['C4']/f**4 + p['C3']/f**3 + p['C2']/f**2
    # Paragraph after Eq. F4, bertaina2024; QKD.ipynb G1 normalization (unrounded).
    g0 = (2*np.pi*p['B_hz'])**2*(1+p['delta'])/(1+p['gamma'])
    # Eq. F4, bertaina2024.
    gain = g0/(2j*np.pi*f)**2*(1j*f+p['B_hz']*p['gamma'])/(1j*f+p['B_hz']*p['delta'])
    # Eq. F2, bertaina2024.
    laser = cav + free/abs(1+gain)**2 if cavity else free
    # Eq. 6, bertaina2024.
    fiber_free = p['l']*(a+b)/f**2*(p['fc1_hz']/(f+p['fc1_hz']))**2
    # Eq. 8, bertaina2024: length-dependent term; detector excluded, Appendix G.
    fiber_dual = ((p['lambda_s_nm']-p['lambda_q_nm'])/p['lambda_s_nm'])**2*p['l']*(a+b)/f**2
    fiber = fiber_dual if stabilized else fiber_free
    # Appendix G, bertaina2024, unnumbered S_detection: one central term, no arm factor.
    det = p['s0']*(p['fc2_hz']/(f+p['fc2_hz']))**2 if stabilized else 0.
    # Eq. 5, bertaina2024; supplementary check for the other Table I rows.
    if common:
        return 4*np.sin(2*np.pi*f*p['n']*(a-b)/p['c_km_s'])**2*laser + p['K']*fiber + det
    # Eq. 7, bertaina2024; primary independent-laser audit (Table I 6/7).
    return 2*laser + fiber + det


def calculate(config, scenario):
    p, op, grid = (config[key] for key in ('physics','operation','grid'))
    # Eq. 5 geometry definition, bertaina2024 Table I: fixed LB, LA=LB+deltaL.
    b = op['LB_km']; a = b + scenario['delta_L_km']
    def psd(f):
        return spectra(f,p,a,b,common=scenario['common'],stabilized=scenario['stabilized'],cavity=scenario['cavity'])
    def integral(lower):
        # Eq. 4, bertaina2024: logarithmic change of integration variable, df=exp(x) dx.
        return quad(lambda x:psd(np.exp(x))*np.exp(x),np.log(lower),np.log(grid['f_max_hz']),
                    epsabs=QUADRATURE_RTOL*op['sigma_limit_rad']**2,
                    epsrel=QUADRATURE_RTOL,limit=QUADRATURE_SUBDIVISIONS)[0]
    # Eq. 4 / Sec. V, bertaina2024: solve threshold; tau_max is a cap, not working time.
    cap_variance = integral(1/op['tau_max_s'])
    if cap_variance <= op['sigma_limit_rad']**2:
        tau = op['tau_max_s']
    else:
        cutoff = brentq(lambda f:integral(f)-op['sigma_limit_rad']**2,1/op['tau_max_s'],grid['f_max_hz'],
                        xtol=ROOT_XTOL_HZ,rtol=ROOT_RTOL)
        tau = 1/cutoff
    # Eq. 4 and unnumbered duty formula, Sec. I p.3, bertaina2024.
    return dict(scenario=scenario['name'],variance_rad2=integral(1/tau),tau_q_s=tau,
                duty=tau/(tau+op['tau_ps_s']),cap_variance_rad2=cap_variance)


def authors(notebook):
    namespace = dict(np=np,math=math)
    nb = json.loads(Path(notebook).read_text())
    # QKD.ipynb Cells 13,15,17,19,21: AST extraction preserves all protocol functions verbatim.
    for cell in (13,15,17,19,21):
        tree = ast.parse(''.join(nb['cells'][cell]['source']))
        tree.body = [node for node in tree.body if isinstance(node,ast.FunctionDef)]
        exec(compile(tree,f'QKD.ipynb:cell{cell}:unmodified','exec'),namespace)
    return namespace


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',default=str(ROOT/'configs/config.toml'))
    parser.add_argument('--output',default='results/audit/independent.json');args=parser.parse_args()
    with Path(args.config).open('rb') as stream:c=tomllib.load(stream)
    rows=[calculate(c,s) for s in c['scenarios']]
    # Primary requested case is row 7; row 6 tests independent lasers without dual compensation.
    Path(args.output).write_text(json.dumps(dict(primary_scenario='7',auxiliary_scenarios=['1','2','3','4','5','6'],
        method='adaptive log-frequency quadrature + bracketed cutoff root; no core import',rows=rows),indent=2))
    print(json.dumps(rows,indent=2))
