"""Part C report rendering, standard-library HTML and matplotlib scientific plots.

Run: python -m tfqkd.lab_report RESULT.json --output DIRECTORY
Numeric results are already computed; rendering does not alter calculations.
"""
from __future__ import annotations
import html
import base64
import json
from pathlib import Path
import numpy as np

def number(value):
    if value is None:return 'undefined'
    if isinstance(value,str):return value
    return f'{value:.9g}'

def table(headers,rows):
    # Markdown escaping is presentation only. HTML is escaped independently when rendered.
    cell=lambda value:str(value).replace('|','\\|').replace('\n',' ')
    return ['| '+' | '.join(map(cell,headers))+' |','| '+' | '.join('---' for _ in headers)+' |']+['| '+' | '.join(cell(value) for value in row)+' |' for row in rows]

def receiver_lines(selected):
    r=selected.get('receiver')
    if not r:return []
    lines=['','### Losses and receiver representation','',r['loss_convention'],
           'Complete arm losses [dB]: '+str(r['complete_arm_losses_db'])+'; effective equalized loss [dB]: '+number(r['effective_equalized_loss_db'])+'.',
           '',r['reduction']]
    if r['channels']:
        lines+=['']+table(['Detector projection','Efficiency','Dark [Hz]','Explicit background [Hz]','Protocol noise [Hz]','Key [bit/s]','Raw secure-rate bound [bit/s]'],
            [(c['name'],number(c['efficiency']),number(c['dark_count_rate_hz']),number(c.get('background_count_rate_hz')),number(c['protocol_noise_count_rate_hz']),number(c['key_bps']),number(c['raw_key_bps'])) for c in r['channels']])
        lines+=['','Scalar projection interval [bit/s]: '+str(r['projection_range_bps'])+'. This interval is not a physical two-detector confidence interval.']
    return lines

def target_lines(result):
    lines = []
    requirements = result['selected'].get('protocol_requirements')
    if requirements:
        lines += ['', '## Conditional protocol requirements', '', requirements['status'],
                  requirements['interpretation'], 'Explicit target [bit/s]: '+number(requirements['target_key_bps'])+'.']
        if requirements['requirements']:
            lines += ['']+table(['Parameter','Detector projection','Maximum','Conditions','Verified key [bit/s]','Status'],
                [(r['parameter'],r['detector_projection'],number(r['maximum']),
                  'f_EC='+number(r['conditional_f_ec']) if r['conditional_f_ec'] is not None else 'e_d='+number(r['conditional_e_d']),
                  number(r.get('verified_key_bps')),r['status']) for r in requirements['requirements']])
            lines += ['', 'The e_d requirement uses f_EC='+number(requirements['reference_f_ec'])+
                      ' as a conditional literature reference, not the main calculation default: '+requirements['reference_source']+'.']
    reach = result.get('reach')
    if reach and reach['status'] == 'disabled':
        return lines + ['', '## Reach and rate uncertainty', '', reach['message']]
    if reach:
        lines += ['', '## Reach and rate uncertainty', '', reach['message'], reach['status'],
                  'Zero-rate reach [total km]: '+number(reach['zero_rate_reach_km'])+'.']
        if reach.get('target_reach'):
            lines += ['']+table(['Case','Threshold [bit/s]','Reach [total km]','Status'],
                [(k,number(v['threshold_bps']),number(v['reach_km']),v['status']) for k,v in reach['target_reach'].items()])
        if reach.get('target_status'):lines += ['',reach['target_status']]
        if reach.get('rate_curve_basis'):lines += ['',reach['rate_curve_basis']]
        if reach.get('endpoint_scheduling'):lines += ['',reach['endpoint_scheduling']]
        if reach.get('limitation'):lines += ['',reach['limitation']]
        lines += ['', 'A rate multiplier leaves the mathematical zero unchanged. Kilometre shifts above refer to the explicit target rate and the computed length curve.']
    return lines

def observed_markdown(result):
    s=result['selected'];phase=result['mode']=='measured_phase'
    lines=[result['message'],'','## 1. Phase input and key estimate','',
           f"Protocol: {s['protocol']}; reference key rate: {number(s['key_bps'])} bit/s; raw bound: {number(s['raw_key_bps'])} bit/s.",
           f"sigma_phi: {number(s['sigma_phi_rad'])} rad; tau_Q: {number(s['tau_q_s'])} s; duty: {number(s['duty'])}.",s['rate_estimate']+'.']
    if phase:
        lines+=['',f"tau_PS: {number(s['phase_measurement']['tau_ps_s'])} s; e_phi (Eq.1, bertaina2024): {number(s['e_phi'])}.",
                'The supplied window is used directly; no phase-threshold solve, spectral integration, stabilization ceiling, compensation comparison, loop recommendation or spectral decomposition is computed.',
                'The protocol retains the unchanged author Gaussian phase-error implementation; displayed e_phi=sigma_phi^2/4 is the small-phase Eq.1 approximation.',
                'Observation / explicit-bound metadata: '+json.dumps(s['phase_measurement'],ensure_ascii=False)]
    else:
        lines+=['','Direct PSD nodes: '+', '.join(result['direct_nodes'])+'.',
                'No inverse coefficient requirements or model-coefficient sensitivity is inferred for those nodes. Select mode=fit to identify a parametric model for that analysis.']
        if result['fiber_ceiling']:
            ceiling=result['fiber_ceiling']
            lines+=['',ceiling['message'],'', '## 2. Fiber stabilization ceiling','']
            lines+=table(['Model','Key [bit/s]'],[(k,number(v)) for k,v in ceiling['rates_bps'].items()])
            lines+=['','Realized fractions: '+json.dumps(ceiling['realized_ceiling_fraction'])+'.']
        else:lines+=['',result['comparison_notice']]
    lines+=receiver_lines(s)
    if not phase:
        analysis=result['sensitivity']
        lines+=['','## Ranked sensitivity for remaining explicit parameters','',analysis['interpretation'],'']
        lines+=table(['Parameter','Applied change','Key [bit/s]','Gain [bit/s]'],[(r['parameter'],r['operation'],number(r['key_bps']),number(r['gain_bps'])) for r in analysis['ranked']])
        lines+=['','Not evaluated:','']+table(['Parameter','Reason'],[(r['parameter'],r['status']) for r in analysis['not_ranked']])
        for row in analysis['ranked']:
            if row.get('note'):lines+=['',row['parameter']+': '+row['note']]
        d=result['variance_contributions'];lines+=['','## Phase-variance contributions','']
        lines+=table(['Contribution','Variance [rad^2]','Fraction'],[(r['component'],number(r['variance_rad2']),number(r['fraction'])) for r in d['rows']])
        lines+=['','![Phase-noise spectrum](spectrum.png)','', '![Accumulated variance](cumulative_variance.png)', '', 'Dense plots retain endpoints and extrema in each log-frequency bin. Full numerical spectra and accumulated integrals are in result.json.']
        if result.get('input_ranges'):
            lines+=['','Range results: '+json.dumps(result['input_ranges']['output_ranges']),result['input_ranges']['range_method']]
    requirements=s.get('equipment_requirements',{}).get('requirements',[])
    if requirements:
        lines+=['','## Conditional equipment requirements','']+table(['Parameter','Maximum','Status'],[(r['parameter'],number(r['maximum']),r['status']) for r in requirements])
    lines+=target_lines(result)
    lines+=['','## Applicability and Tool limitations','',result['tool_limitations'],
            '', '## Input provenance','']
    lines+=table(['Parameter','Origin'],[(r['parameter'],json.dumps(r,ensure_ascii=False)) for r in s['provenance']])
    lines+=['','## Run notices','']+['- '+w for w in dict.fromkeys(s['warnings'])]
    lines+=['',f"Executed {result['input_case_calculations']} input cases and {result['spectral_calculations']} spectral calculations in {number(result['elapsed_seconds'])} s."]
    return '\n'.join(lines)+'\n'

def _model_markdown(result):
    if result.get('mode')=='direct_psd' and result.get('fiber_ceiling'):
        # A measured laser with a parametric arm model retains the full C report order and loop diagnostics.
        display=dict(result);display.pop('mode')
        display['message']=result['fiber_ceiling']['message']+'\n\n'+result['message']+' Direct nodes: '+', '.join(result['direct_nodes'])+'.'
        if result.get('upper_estimate_notice'):
            display['message']=result['upper_estimate_notice']+'\n\n'+display['message']
        return _model_markdown(display)
    if result.get('mode') in ('measured_phase','direct_psd'):return observed_markdown(result)
    s=result['selected'];ceiling=result['fiber_ceiling'];analysis=result['sensitivity'];decomposition=result['variance_contributions']
    lines=[result['message'],'','## 1. Fiber stabilization ceiling','']
    lines+=table(['Model','Key rate [bit/s]','Estimate'],[(name,number(rate),'optimistic: phase-detection noise omitted' if name=='R_classical' else 'conditional on explicit/optimistic input') for name,rate in ceiling['rates_bps'].items()])
    lines+=['',f"H = {number(ceiling['H'])}; realized classical fraction = {number(ceiling['realized_ceiling_fraction']['classical'])} (optimistic); dual fraction = {number(ceiling['realized_ceiling_fraction']['dual'])}.",
            ceiling['detection_convention']+'.','']
    lines+=['## 2. Selected installation','',f"{s['compensation']}, {s['protocol']}; arms {number(s['arm_lengths_km'][0])} / {number(s['arm_lengths_km'][1])} km.",
            f"Key rate: {number(s['key_bps'])} bit/s; tau_Q: {number(s['tau_q_s'])} s; duty: {number(s['duty'])}.",
            'Estimate: '+s['rate_estimate']+'.']
    if result['input_ranges']:
        ranges=result['input_ranges']
        lines+=['', 'Range output:', '']+table(['Output','Minimum','Maximum'],[(key,number(lo),number(hi)) for key,(lo,hi) in ranges['output_ranges'].items()])
        lines+=['',ranges['range_method']+'.','Headline H and sensitivity use the explicit point inputs; each evaluated range coordinate retains its four-rate result in JSON.']
    lines+=receiver_lines(s)
    lines+=target_lines(result)
    lines+=['','## 3. Ranked sensitivity','',analysis['interpretation'],'']
    lines+=table(['Rank','Parameter','Applied change','New rate [bit/s]','Gain [bit/s]','Relative gain'],[
        (row['rank'],row['label']+' ['+row['units']+']',f"{number(row['display_baseline_value'])} → {number(row['display_changed_value'])} ({row['operation']})",number(row['key_bps']),number(row['gain_bps']),number(row['gain_fraction'])) for row in analysis['ranked']])
    for row in analysis['ranked']:
        if row.get('note'):lines+=['',row['parameter']+': '+row['note']]
    lines+=['','Not ranked:','']+table(['Parameter','Reason'],[(row['parameter'],row['status']+(': '+row['reason'] if 'reason' in row else '')) for row in analysis['not_ranked']])
    lines+=['','Unmeasured amplitudes have no calculated sensitivity; their conditional equipment requirements are below.',
            '', '## 4. Conditional equipment requirements','']
    requirements=s.get('equipment_requirements',{}).get('requirements',[])
    if requirements:
        from .lab_analysis import PARAMETER_LABELS
        def requirement_label(path):
            label,units=PARAMETER_LABELS.get(path,(path,'model units'))
            return label+' ['+units+']'
        lines+=table(['Parameter','Maximum','Verified key loss [fraction]','Status'],[(requirement_label(row['parameter']),number(row['maximum']),number(row.get('verified_loss_fraction')),row['status']) for row in requirements])
        lines+=['','Each limit assumes the other unmeasured contributions are zero; the individual maxima are not joint noise-budget limits.']
    else:lines+=['No missing active noise amplitudes; no inverse requirement is substituted for measured input.']
    lines+=['','## 5. Recommended classical loop','']
    loop=ceiling['results']['classical']['loop']
    lines+=table(['g recommended [s^-1]','g_crit [s^-1]','1-g/g_crit','Safety fraction'],[(number(loop['g_per_s']),number(loop['g_crit_per_s']),number(loop['stability_margin_fraction']),number(loop['safety_fraction']))])
    lines+=['','Recommendation refers to the optimized classical comparison; it does not certify hardware stability. The margin is a gain fraction, not a phase margin in degrees.']
    if s['compensation']=='classical' and s['loop']['mode']=='manual':
        lines+=['',f"Selected manual g = {number(s['loop']['g_per_s'])}; selected margin = {number(s['loop']['stability_margin_fraction'])}. The recommendation above is computed separately."]
    lines+=['','## 6. Phase-variance contributions','',
            f"At the same tau_Q: lower Fourier limit {number(decomposition['lower_hz'])} Hz; total variance {number(decomposition['total_variance_rad2'])} rad^2.",'']
    lines+=table(['Component','Variance [rad^2]','Share'],[(row['component'],number(row['variance_rad2']),number(row['fraction'])) for row in decomposition['rows']])
    lines+=['',decomposition['detection_status']]
    lines+=['',f"Dominant component: {decomposition['dominant_component']}; contribution-sum relative error: {number(decomposition['sum_relative_error'])}.",
            '', '![Phase-noise spectrum](spectrum.png)','', '![Accumulated variance](cumulative_variance.png)', '', 'Dense plots retain endpoints and extrema in each log-frequency bin. Full numerical spectra and accumulated integrals are in result.json.',
            '', '## 7. Applicability and Tool limitations','',result['tool_limitations'],
            '', 'Input provenance (explicit inputs or accepted fits; no equipment presets):','']
    lines+=table(['Parameter','Value / origin'],[(row['parameter'],json.dumps(row,ensure_ascii=False)) for row in s['provenance']])
    diagnostic=ceiling['results']['classical']['classical_diagnostic']
    lines+=['','Secondary bandwidth diagnostic: '+json.dumps(diagnostic,ensure_ascii=False),
            '', 'Run notices:','']
    warnings=s['warnings']+ceiling['warnings']+(result['input_ranges']['warnings'] if result['input_ranges'] else [])
    lines+=['- '+warning for warning in dict.fromkeys(warnings)]
    lines+=['',f"Executed {result['input_case_calculations']} input cases and {result['spectral_calculations']} spectral candidate calculations in {number(result['elapsed_seconds'])} s; report rendering is timed separately when written."]
    return '\n'.join(lines)+'\n'

def markdown(result):
    text = _model_markdown(result)
    web = result.get('web_input')
    if not web:return text
    lines=['', '## Web input conventions', '',
        'SNS intensities: per user, as usually published; core uses twice this value (Appendix D, bertaina2024). '
        'CAL intensity is already per user and is not doubled (Appendix B / author Cell 19).', '']
    lines+=table(['Parameter','Per user [photons/pulse]','Core value','Factor','Source'],
                 [(r['parameter'],number(r['per_user']),number(r['core']),r['factor'],r['source']) for r in web['intensities']])
    if web.get('source'):lines+=['',str(web['source'])]
    if web['ignored_attenuation']:lines+=['','Complete arm losses override '+', '.join(web['ignored_attenuation'])+'; per-km attenuation was not used.']
    if web['warnings']:lines+=['','Input warnings:','']+['- '+path+': '+notice for path,notice in web['warnings'].items()]
    return text+'\n'.join(lines)+'\n'


def html_document(markdown_text,*,title='Twin-Field QKD laboratory report',images=None):
    """Render our emitted Markdown subset; arbitrary input is escaped, no third-party parser."""
    images=images or {};output=[];lines=markdown_text.splitlines();index=0
    while index<len(lines):
        line=lines[index]
        if line.startswith('| ') and index+1<len(lines) and lines[index+1].startswith('| ---'):
            records=[line];index+=2
            while index<len(lines) and lines[index].startswith('| '):records.append(lines[index]);index+=1
            output.append('<div class="table-wrap"><table>')
            for n,row in enumerate(records):
                tag='th' if n==0 else 'td'
                cells=row.strip('| ').replace('\\|','∣').split('|')
                output.append('<tr>'+''.join(f'<{tag}>'+html.escape(cell.strip())+f'</{tag}>' for cell in cells)+'</tr>')
            output.append('</table></div>');continue
        if line.startswith('![') and '](' in line:
            alt,path=line[2:].split('](',1);path=path.rstrip(')')
            output.append(images.get(path,'<img alt="'+html.escape(alt,quote=True)+'" src="'+html.escape(path,quote=True)+'">'))
        elif line.startswith('## '):output.append('<h2>'+html.escape(line[3:])+'</h2>')
        elif line.startswith('# '):output.append('<h1>'+html.escape(line[2:])+'</h1>')
        elif line.startswith('- '):output.append('<p class="notice">'+html.escape(line[2:])+'</p>')
        elif line:output.append('<p>'+html.escape(line)+'</p>')
        index+=1
    # Report styling only; no physical or computational constants.
    style='body{font:16px/1.5 system-ui,sans-serif;max-width:1100px;margin:32px auto;padding:0 20px;color:#182534;background:#fff}h2{margin-top:40px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px 12px;text-align:left;border-bottom:1px solid #ccd6e0;vertical-align:top;overflow-wrap:anywhere}th{background:#edf3f7}.table-wrap{overflow-x:auto}p{overflow-wrap:anywhere}svg,img{max-width:100%;height:auto}.notice{padding-left:12px;border-left:3px solid #899bab}@media print{body{margin:0;max-width:none}h2{break-after:avoid}tr{break-inside:avoid}}'
    return '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>'+html.escape(title)+'</title><style>'+style+'</style><main><h1>'+html.escape(title)+'</h1>'+''.join(output)+'</main></html>\n'

# Presentation-only density; numerical spectra and integrals are kept in full in result.json.
DISPLAY_LOG_BINS = 1024

def display_indices(frequency, values, bins=DISPLAY_LOG_BINS):
    if len(frequency) <= 4*bins:return np.arange(len(frequency))
    edges=np.searchsorted(frequency,np.geomspace(frequency[0],frequency[-1],bins+1))
    indices=[0,len(frequency)-1]
    for left,right in zip(edges[:-1],edges[1:]):
        if right<=left:continue
        data=values[left:right]
        indices.extend((left,right-1,left+int(np.argmin(data)),left+int(np.argmax(data))))
    return np.unique(indices)

def plots(result,directory):
    import os
    from .config import ROOT
    cache=ROOT/'results/.matplotlib'
    cache.mkdir(parents=True,exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR',str(cache))
    import matplotlib
    # Reports export files; GUI backends are unsafe in the local HTTP worker thread.
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    d=result['variance_contributions'];f=np.array(d['frequency_hz']);images={}
    for filename,quantities,ylabel in [('spectrum',d['phase_psd_rad2_per_hz'],'One-sided phase PSD [rad²/Hz]'),('cumulative_variance',d['cumulative_variance_rad2'],'Variance above lower Fourier limit [rad²]')]:
        fig,ax=plt.subplots(figsize=(8,4.5))
        for name,values in quantities.items():
            y=np.asarray(values);valid=y>0
            if np.any(valid):
                indices=display_indices(f[valid],y[valid])
                ax.loglog(f[valid][indices],y[valid][indices],label=name)
        ax.set_xscale('log');ax.set_xlim(f[0],f[-1])
        if not any(np.any(np.asarray(values)>0) for values in quantities.values()):
            ax.text(.5,.5,'All contributions are zero under these inputs',ha='center',transform=ax.transAxes)
        ax.axvline(d['lower_hz'],color='black',linestyle=':',label='1 / tau_Q')
        ax.set(xlabel='Fourier frequency [Hz]',ylabel=ylabel);ax.legend();ax.grid(alpha=.3);fig.tight_layout()
        fig.savefig(directory/(filename+'.png'))
        plt.close(fig)
        encoded=base64.b64encode((directory/(filename+'.png')).read_bytes()).decode('ascii')
        images[filename+'.png']='<img alt="'+html.escape(filename.replace('_',' '),quote=True)+'" src="data:image/png;base64,'+encoded+'">'
    return images

def write_report(result,directory):
    from time import perf_counter
    started=perf_counter();directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    images={} if result.get('mode')=='measured_phase' else plots(result,directory)
    text=markdown(result)
    (directory/'report.md').write_text(text)
    (directory/'report.html').write_text(html_document(text,images=images))
    result['plot_display'] = 'For dense grids: plotted samples retain bin endpoints and minimum/maximum in each log-frequency bin; full numerical series remains in JSON. This is display reduction only.'
    result['report_render_seconds']=perf_counter()-started
    (directory/'result.json').write_text(json.dumps(result,separators=(',',':'),allow_nan=False))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('result',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();write_report(json.loads(args.result.read_text()),args.output)
