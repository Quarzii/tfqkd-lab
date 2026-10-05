"""Compare two/three explicitly specified setups and their working-length outputs.

Working lengths mean physical TOTAL length; each variant keeps its own signed
imbalance. Default coordinates come from the user's C2 specification, not a preset.
Run: python -m tfqkd.lab_compare COMPARISON.toml --output DIRECTORY
"""
from __future__ import annotations
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from .config import ROOT
from .lab_inputs import resolve,read_toml,validate
from .lab_engine import calculate_point
from .engine import independent_jobs
from .lab_run import run
from .lab_report import number,table,html_document,write_report

def _variant(job):
    label,path=job
    # Independent setups are parallelized here; avoid nested processes in each sensitivity/range evaluation.
    return dict(label=label,configuration=str(path),result=run(path,workers=1))

def _working(job):
    label,lab,total=job
    a,b=lab.arms()
    # Eq.5 geometry, bertaina2024: physical total scan keeps the variant's prescribed imbalance.
    lab.settings['line']=dict(length_km=total,imbalance_km=a-b);validate(lab)
    nominal=calculate_point(lab);ranges=None
    if lab.ranges:
        # Arm-length/total ranges cannot describe a fixed nominal working-length scan without a new geometry rule.
        geometry=[key for key in lab.ranges if key.startswith('line.') and key!='line.imbalance_km']
        if geometry:raise ValueError('Working-length comparison requires a fixed total: remove total/arm-length ranges or specify separate variants')
        from .lab_uncertainty import calculate_ranges
        ranges=calculate_ranges(lab,workers=1,calculator=calculate_point)
    return dict(label=label,physical_total_km=total,arm_lengths_km=list(lab.arms()),key_bps=nominal['key_bps'],
                tau_q_s=nominal['tau_q_s'],duty=nominal['duty'],variance_rad2=nominal['variance_rad2'],
                compensation=nominal['compensation'],g_per_s=nominal['loop']['g_per_s'],g_crit_per_s=nominal['loop']['g_crit_per_s'],
                upper_estimate=bool(lab.missing_amplitudes) or bool(lab.settings.get('_ideal_protocol_fields')) or nominal['compensation']=='classical',
                output_ranges=ranges['output_ranges'] if ranges else None,
                input_case_calculations=1+(ranges['input_case_calculations'] if ranges else 0),
                spectral_calculations=nominal['spectral_calculations']+(ranges['spectral_calculations'] if ranges else 0))

def compare(user,*,workers=None):
    started=perf_counter()
    raw=read_toml(user) if isinstance(user,(str,Path)) else user
    base=Path(user).resolve().parent if isinstance(user,(str,Path)) else Path.cwd()
    if set(raw)!={'comparison'}:raise ValueError('Comparison file must contain only [comparison]; setups are fully specified separate files')
    spec=raw['comparison'];allowed={'configurations','labels','working_total_lengths_km','workers'}
    if set(spec)-allowed:raise ValueError('Unknown comparison fields: '+', '.join(sorted(set(spec)-allowed)))
    configurations=spec.get('configurations',[]);labels=spec.get('labels',[])
    if len(configurations) not in (2,3) or len(labels)!=len(configurations) or len(set(labels))!=len(labels):
        raise ValueError('Supply two or three configurations with distinct labels of matching length')
    if not all(isinstance(label,str) and label.strip() for label in labels):raise ValueError('Comparison labels must be nonempty strings')
    paths=[(base/path).resolve() for path in configurations]
    labs=[resolve(path) for path in paths]
    from .lab_observed import measured_phase,explicit_losses
    for label,lab in zip(labels,labs):
        if measured_phase(lab):
            raise ValueError(f'{label}: measured phase is a point observation; compensation and working-length forecasts are unavailable')
        if 'line' in lab.direct_spectra:
            raise ValueError(f'{label}: a direct arm residual does not identify other compensation states or working-length forecasts')
        if explicit_losses(lab.settings['line']) is not None:
            raise ValueError(f'{label}: complete losses are point measurements; working-length comparison requires independently supplied losses at each length, not an inferred attenuation coefficient')
    count=spec.get('workers',labs[0].settings['numerics']['workers']) if workers is None else workers
    if not isinstance(count,int) or count<1:raise ValueError('Comparison workers must be a positive integer')
    defaults=read_toml(ROOT/'configs/lab_defaults.toml')['comparison']
    lengths=spec.get('working_total_lengths_km',defaults['working_total_lengths_km'])
    if not isinstance(lengths,list) or not lengths or not np.all(np.isfinite(lengths)) or min(lengths)<=0 or len(set(lengths))!=len(lengths):
        raise ValueError('Working total lengths must be a nonempty list of distinct positive km values')
    variants=independent_jobs(_variant,list(zip(labels,paths)),count)
    jobs=[(label,lab.clone(),float(length)) for length in lengths for label,lab in zip(labels,labs)]
    rows=independent_jobs(_working,jobs,count)
    for length in lengths:
        at_length=[row for row in rows if row['physical_total_km']==length];baseline=at_length[0]['key_bps']
        for row in at_length:
            # User-requested working-rate ratio; zero denominator cases remain explicit.
            row['rate_ratio_to_first']=row['key_bps']/baseline if baseline>0 else ('infinite' if row['key_bps']>0 else None)
            row['tau_ratio_to_first']=row['tau_q_s']/at_length[0]['tau_q_s']
            row['duty_ratio_to_first']=row['duty']/at_length[0]['duty'] if at_length[0]['duty']>0 else None
    return dict(mode='comparison',variants=variants,working_lengths=rows,working_total_lengths_km=lengths,
                ratio_baseline_label=labels[0],elapsed_seconds=perf_counter()-started,
                input_case_calculations=sum(v['result']['input_case_calculations'] for v in variants)+sum(r['input_case_calculations'] for r in rows),
                spectral_calculations=sum(v['result']['spectral_calculations'] for v in variants)+sum(r['spectral_calculations'] for r in rows),
                interpretation='Explicit setups; no preset values. Working lengths are total physical lengths, with each signed imbalance fixed. Range comparisons retain sampled-envelope limitations.')

def comparison_markdown(result):
    from .lab_limits import CLASSICAL_DETECTION_NOISE,tool_limitations_markdown
    lines=['# Installation comparison','',result['interpretation'],'', '## Configured points','']
    ideal = [v['label'] for v in result['variants'] if v['result'].get('upper_estimate_notice')]
    if ideal:
        lines = ['Upper model estimates: unreported e_d/f_EC evaluated at ideal limits in '+', '.join(ideal)+'.', '']+lines
    lines+=table(['Variant','Total [km]','Scheme','Key [bit/s]','tau_Q [s]','Duty','H','Classical fraction (optimistic)','Dual fraction'],[
        (v['label'],number(v['result']['selected']['physical_length_km']),v['result']['selected']['compensation'],
         number(v['result']['selected']['key_bps']),number(v['result']['selected']['tau_q_s']),number(v['result']['selected']['duty']),
         number(v['result']['fiber_ceiling']['H']),number(v['result']['fiber_ceiling']['realized_ceiling_fraction']['classical']),
         number(v['result']['fiber_ceiling']['realized_ceiling_fraction']['dual'])) for v in result['variants']])
    lines+=['','## Sensitivity and loop recommendations','']
    summaries=[]
    for variant in result['variants']:
        r=variant['result'];ranked=r['sensitivity']['ranked'];winner=ranked[0] if ranked else None
        loop=r['fiber_ceiling']['results']['classical']['loop'];dominant=r['variance_contributions']['dominant_component']
        summaries.append((variant['label'],winner['label'] if winner else 'none',number(winner['gain_fraction'] if winner else None),
                          number(loop['g_per_s']),number(loop['g_crit_per_s']),number(loop['stability_margin_fraction']),dominant,
                          ', '.join(r['selected']['missing_amplitudes']) or 'none'))
    lines+=table(['Variant','Largest tested gain','Relative gain','g recommended [s^-1]','g_crit [s^-1]','Margin','Dominant noise','Unmeasured: see requirements'],summaries)
    lines+=['','## Working total lengths','',f"Rate ratios use {result['ratio_baseline_label']} as the denominator.",'']
    lines+=table(['Total [km]','Variant','Scheme','Key [bit/s]','tau_Q [s]','Duty','Rate ratio','tau ratio','Duty ratio','Estimate'],[
        (number(row['physical_total_km']),row['label'],row['compensation'],number(row['key_bps']),number(row['tau_q_s']),number(row['duty']),number(row['rate_ratio_to_first']),number(row['tau_ratio_to_first']),number(row['duty_ratio_to_first']),
         'optimistic' if row['upper_estimate'] else 'explicit-input model point') for row in result['working_lengths']])
    lines+=['',CLASSICAL_DETECTION_NOISE,'','Detailed ranked sensitivity, inverse requirements, loop recommendation, contribution plots and provenance are in each variant report.']
    for i,v in enumerate(result['variants'],start=1):lines+=['',f"Variant {i}: {v['label']} — report in variant_{i}/report.html and variant_{i}/report.md."]
    lines+=['','## Tool limitations','',tool_limitations_markdown(),'',f"Executed {result['input_case_calculations']} input cases, {result['spectral_calculations']} spectral candidates; compute time {number(result['elapsed_seconds'])} s."]
    return '\n'.join(lines)+'\n'

def write_comparison(result,directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    for i,variant in enumerate(result['variants'],start=1):write_report(variant['result'],directory/f'variant_{i}')
    text=comparison_markdown(result)
    (directory/'report.md').write_text(text);(directory/'report.html').write_text(html_document(text,title='Twin-Field QKD installation comparison'))
    (directory/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=compare(args.config);write_comparison(result,args.output);print(comparison_markdown(result))
