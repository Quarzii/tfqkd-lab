"""Reproducible B6 digitization and diagnostic fits, authorized by the user.

Poppler renders are retained with calibrated log-axis anchors and exact extracted
pixel coordinates. Log least squares is an adopted numerical method, not an
equation or parameter-estimation method attributed to Bertaina. Conditional fit
SE describes point scatter; it is NOT a device/population uncertainty.
Run: python -m tfqkd.published_fit
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .config import ROOT
from .spectra import free_fiber, free_laser

OUT=ROOT/'results/stage3_b_amendments/b6'

def digitize(spec, settings):
    from matplotlib.image import imread
    im=imread(ROOT/spec['image'])[:,:,:3]
    left,right=spec['x_pixels'];top,bottom=spec['y_pixels']
    red,green,blue=np.moveaxis(im,-1,0)
    if spec['color']=='red':
        mask=(red-green>settings['color_difference'])&(red-blue>settings['color_difference'])
    elif spec['color']=='blue':
        mask=(blue-red>settings['color_difference'])&(blue-green>settings['blue_green_difference'])
    else:
        mask=(np.max(im,axis=-1)<settings['black_threshold'])&(np.ptp(im,axis=-1)<settings['gray_difference'])
    for x0,x1,y0,y1 in spec.get('exclude_rectangles',[]):mask[y0:y1+1,x0:x1+1]=False
    pixels=[]
    for x in range(left+settings['border_pixels'],right-settings['border_pixels'],settings['pixel_step']):
        ys=np.flatnonzero(mask[top+settings['border_pixels']:bottom-settings['border_pixels'],x])+top+settings['border_pixels']
        if spec.get('exclude_white_linewidth_reference'):
            # Fig. 1/Fig. 3, jiang2008: dashed 1/(pi*f^2) reference identifies the ONE-SIDED normalization.
            logf=spec['x_log10'][0]+(x-left)/(right-left)*np.ptp(spec['x_log10'])
            logp=-np.log10(np.pi)-2*logf
            reference_y=top+(spec['y_log10'][1]-logp)/np.ptp(spec['y_log10'])*(bottom-top)
            ys=ys[abs(ys-reference_y)>settings['reference_exclusion_pixels']]
        if not len(ys):continue
        # Digitization only: cluster antialiased line pixels; upper black trace is selected in Jiang Fig. 3.
        clusters=np.split(ys,np.flatnonzero(np.diff(ys)>1)+1)
        selected=clusters[0] if spec.get('upper_black_trace') else ys
        pixels.append((x,float(np.median(selected))))
    xy=np.asarray(pixels)
    # Published figure axis calibration, B6 adopted log-axis digitization; not a physical equation.
    logf=spec['x_log10'][0]+(xy[:,0]-left)/(right-left)*np.ptp(spec['x_log10'])
    logp=spec['y_log10'][1]-(xy[:,1]-top)/(bottom-top)*np.ptp(spec['y_log10'])
    f,p=10**logf,10**logp
    # Eqs. 1,5, didomenico2010; Eq. F1, bertaina2024: single-sided frequency PSD to phase PSD.
    phase=p/f**2 if spec['quantity']=='frequency' else p
    # No two-sided, SSB or round-trip conversion without established metadata. All fitted curves here are one-sided, one-way.
    phase*=spec['normalization_to_one_sided_one_way']
    return xy,f,phase,p

def fit_points(f,phase,spec,settings):
    fiber=spec['model']=='fiber'
    names=('l','fc1_hz') if fiber else ('r3','r2','fc_hz')
    def model(log_parameters):
        values=dict(zip(names,np.exp(log_parameters)))
        # Eq. 6 / Eq. F1, bertaina2024: fit the UNCHANGED source model, no added floor, peak or term.
        return free_fiber(f,spec['length_km'],values) if fiber else free_laser(f,values)
    def residual(parameters):return np.log(model(parameters))-np.log(phase)
    # B6 adopted equal log-frequency/log-residual fit; bounds are floating-point search limits, NOT physical presets.
    if fiber:
        start=np.log([np.median(phase*f*f)/spec['length_km'],np.median(f)])
    else:
        start=np.log([np.median(phase*f**3),np.median(phase*f*f),np.median(f)])
    margin=settings['log_parameter_search_margin']
    lower=np.full(len(names),-margin);upper=np.full(len(names),margin)
    solutions=[least_squares(residual,np.clip(start+np.r_[np.zeros(len(names)-1),shift],lower,upper),
                             bounds=(lower,upper),max_nfev=settings['max_nfev'])
               for shift in settings['cutoff_start_shifts']]
    opt=min(solutions,key=lambda value:value.cost)
    if not opt.success:raise ArithmeticError(opt.message)
    prediction=model(opt.x);ratios=prediction/phase
    _,singular,vh=np.linalg.svd(opt.jac,full_matrices=False)
    # B6 adopted local least-squares SE: residual scatter times inverse information, conditional iid log-point model.
    # This is a numerical fit diagnostic, not measured statistical error bars or source-provided uncertainty.
    formal_rank=np.linalg.matrix_rank(opt.jac)
    # Numerical resolution of the finite-difference Jacobian, derived from machine epsilon (not a physical fit gate).
    rank_tolerance=np.sqrt(np.finfo(float).eps)*max(opt.jac.shape)*singular[0]
    rank=int(np.sum(singular>rank_tolerance))
    # Use the Jacobian SVD directly: squaring its condition number can silently erase a weak direction.
    inverse=np.divide(1.,singular,out=np.zeros_like(singular),where=singular>0)
    covariance=(vh.T*inverse**2)@vh*(2*opt.cost/(len(f)-len(names)))
    parameters=np.exp(opt.x);se=parameters*np.sqrt(np.diag(covariance))
    errors=se.tolist()
    if rank<len(names):
        for idx in np.flatnonzero(np.any(abs(vh[rank:])>np.sqrt(np.finfo(float).eps),axis=0)):errors[idx]=None
    result=dict(parameters=dict(zip(names,parameters.tolist())),conditional_fit_se=dict(zip(names,errors)),
                covariance_log_parameters=covariance.tolist(),jacobian_rank=int(rank),formal_jacobian_rank=int(formal_rank),
                rank_tolerance=float(rank_tolerance),formal_fit_se=dict(zip(names,se.tolist())),singular_values=singular.tolist(),
                rms_log10_residual=float(np.sqrt(np.mean(np.log10(ratios)**2))),
                maximum_multiplicative_discrepancy=float(np.max(np.maximum(ratios,1/ratios))),
                frequency_range_hz=[float(f.min()),float(f.max())],points=len(f),success=bool(opt.success),
                uncertainty_label='conditional local fit SE only; spectral features, digitization correlations and apparatus uncertainty are not represented')
    # B6 identifiability diagnostic: profile only the cutoff over observed span and computational extrapolation factors.
    profile=[]
    for multiplier in settings['cutoff_profile_multipliers']:
        cutoff=f.max()*multiplier
        fixed=np.log(cutoff)
        candidate=least_squares(lambda q:residual(np.r_[q,fixed]),opt.x[:-1],bounds=(lower[:-1],upper[:-1]),max_nfev=settings['max_nfev'])
        profile.append(dict(cutoff_hz=float(cutoff),rms_log10_residual=float(np.sqrt(np.mean(candidate.fun**2))/np.log(10))))
    result['cutoff_profile']=profile
    return result,prediction

def run():
    from matplotlib import pyplot as plt
    manifest=json.loads((ROOT/'sources/data/b6/digitization.json').read_text());settings=manifest['numerics']
    OUT.mkdir(parents=True,exist_ok=True);results=[]
    for spec in manifest['curves']:
        xy,f,phase,plotted=digitize(spec,settings)
        np.savetxt(OUT/(spec['id']+'_points.csv'),np.column_stack((xy,f,plotted,phase)),delimiter=',',
                   header='pixel_x,pixel_y,frequency_hz,plotted_one_sided_psd,one_way_phase_psd_rad2_per_hz',comments='')
        fit,predicted=fit_points(f,phase,spec,settings)
        fig,ax=plt.subplots();ax.loglog(f,phase,'.',label='Digitized source');ax.loglog(f,predicted,label='Eq. 6' if spec['model']=='fiber' else 'Eq. F1')
        ax.set(xlabel='Fourier frequency [Hz]',ylabel='One-sided phase PSD [rad²/Hz]',title=spec['id']);ax.legend();fig.tight_layout();fig.savefig(OUT/(spec['id']+'_fit.png'));plt.close(fig)
        results.append(dict(id=spec['id'],specification=spec,fit=fit))
    (OUT/'fits.json').write_text(json.dumps(dict(method=__doc__,results=results),indent=2))
    return results

if __name__=='__main__':
    for row in run():print(row['id'],row['fit']['parameters'],row['fit']['rms_log10_residual'])
