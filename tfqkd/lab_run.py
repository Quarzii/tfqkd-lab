"""Universal laboratory command line. Inputs are explicit; no preset library.

Run: python -m tfqkd.lab_run INPUT.toml --output DIRECTORY
Part C adds ranked sensitivity, variance contributions and Markdown/HTML output.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from time import perf_counter
from .lab_inputs import resolve, MissingInputs
from .measured_inputs import MeasurementRejected
from .lab_engine import calculate
from .lab_ceiling import stabilization_ceiling,comparison_ready
from .lab_limits import tool_limitations_markdown,limitation_ids

def run_point(lab):
    ceiling=stabilization_ceiling(lab)
    scheme=lab.settings['scheme']['compensation']
    same_point=not (scheme=='classical' and 'g_per_s' in lab.settings['loop'])
    selected=calculate(lab, point=ceiling['results'][{'none':'free','classical':'classical','dual':'dual'}[scheme]] if same_point else None)
    return dict(selected, message=ceiling['message'],fiber_ceiling=ceiling,
                input_case_calculations=selected['input_case_calculations']+ceiling['input_case_calculations'],
                spectral_calculations=selected['spectral_calculations']+ceiling['spectral_calculations'])

def _run(user, *, workers=None):
    started=perf_counter();lab=user if hasattr(user,'settings') else resolve(user)
    from .lab_observed import measured_phase
    if measured_phase(lab):
        point=calculate(lab)
        return dict(mode='measured_phase',message='Estimate from a supplied residual phase observation or explicit phase bound, not a spectral forecast.',
                    selected=point,elapsed_seconds=perf_counter()-started,
                    input_case_calculations=point['input_case_calculations'],spectral_calculations=0,
                    tool_limitations=tool_limitations_markdown(),limits=limitation_ids())
    if lab.direct_spectra:
        from .lab_analysis import sensitivity,variance_contributions
        point=calculate(lab);analysis=sensitivity(lab,point,workers=workers);decomposition=variance_contributions(lab,point)
        ranges=None
        if lab.ranges:
            from .lab_uncertainty import calculate_ranges
            ranges=calculate_ranges(lab,workers=workers,calculator=calculate)
        ceiling=None;notice=None
        if 'line' in lab.direct_spectra:
            notice='No fiber ceiling or compensation-scheme comparison: the uploaded arm PSD describes the selected state and cannot identify a different compensation state or controller.'
        else:
            try:ceiling=stabilization_ceiling(lab)
            except MissingInputs as error:notice='Compensation ceiling omitted: '+str(error)
        return dict(mode='direct_psd',message='Calculation from measured spectral samples; parametric shapes are not imposed on direct nodes.',
                    selected=point,fiber_ceiling=ceiling,comparison_notice=notice,sensitivity=analysis,input_ranges=ranges,
                    variance_contributions=decomposition,direct_nodes=list(lab.direct_spectra),
                    elapsed_seconds=perf_counter()-started,input_case_calculations=point['input_case_calculations']+analysis['input_case_calculations']+(ceiling['input_case_calculations'] if ceiling else 0)+(ranges['input_case_calculations'] if ranges else 0),
                    spectral_calculations=point['spectral_calculations']+analysis['spectral_calculations']+decomposition['spectral_calculations']+(ceiling['spectral_calculations'] if ceiling else 0)+(ranges['spectral_calculations'] if ranges else 0),
                    tool_limitations=tool_limitations_markdown(),limits=limitation_ids())
    lab=comparison_ready(lab)
    point=run_point(lab);ceiling=point.pop('fiber_ceiling');point.pop('message')
    ranges=None
    if lab.ranges:
        from .lab_uncertainty import calculate_ranges
        ranges=calculate_ranges(lab,workers=workers,calculator=run_point)
        # Count physical point cases including four comparisons and inverse evaluations, not just endpoint coordinates.
        ranges['input_case_calculations']=sum(row['result']['input_case_calculations'] for row in ranges['cases'])
    from .lab_analysis import sensitivity,variance_contributions
    analysis=sensitivity(lab,point,workers=workers)
    decomposition=variance_contributions(lab,point)
    return dict(message=ceiling['message'],selected=point,fiber_ceiling=ceiling,input_ranges=ranges,
                sensitivity=analysis,variance_contributions=decomposition,
                elapsed_seconds=perf_counter()-started,
                input_case_calculations=point['input_case_calculations']+(ranges['input_case_calculations'] if ranges else 0)+analysis['input_case_calculations'],
                spectral_calculations=point['spectral_calculations']+(ranges['spectral_calculations'] if ranges else 0)+analysis['spectral_calculations']+decomposition['spectral_calculations'],
                tool_limitations=tool_limitations_markdown(),limits=limitation_ids())

def run(user, *, workers=None):
    started = perf_counter()
    lab = resolve(user)
    result = _run(lab, workers=workers)
    from .lab_reach import reach_report
    ranges = result.get('input_ranges')
    reach = reach_report(lab, range_cases=ranges['cases'] if ranges else None, workers=workers,
                         monotone_parameters=ranges['monotone_parameters'] if ranges else ())
    result['reach'] = reach
    result['input_case_calculations'] += reach['input_case_calculations']
    result['spectral_calculations'] += reach['spectral_calculations']
    if lab.settings.get('_ideal_protocol_fields'):
        missing = ', '.join(lab.settings['_ideal_protocol_fields'])
        result['upper_estimate_notice'] = 'Upper model estimate: unreported '+missing+' evaluated at ideal limiting values; not measured apparatus defaults.'
        result['message'] = result['upper_estimate_notice']+'\n\n'+result['message']
    result['elapsed_seconds'] = perf_counter()-started
    return result

def markdown(result):
    from .lab_report import markdown as render
    return render(result)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('config',type=Path);parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    try:
        from .lab_inputs import read_toml
        if 'comparison' in read_toml(args.config):
            from .lab_compare import compare,write_comparison
            result=compare(args.config)
            if args.output:write_comparison(result,args.output)
            from .lab_compare import comparison_markdown
            print(comparison_markdown(result));return 0
        result=run(args.config)
    except (MissingInputs,MeasurementRejected,ValueError,ArithmeticError,OSError) as error:
        print(json.dumps(dict(error=str(error),missing_inputs=getattr(error,'parameters',None),diagnostics=getattr(error,'diagnostics',None)),indent=2));return 2
    if args.output:
        from .lab_report import write_report
        write_report(result,args.output)
    print(markdown(result));return 0

if __name__=='__main__':raise SystemExit(main())
