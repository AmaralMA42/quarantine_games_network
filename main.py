"""Single entry point for Quarantine Games dynamics, sweeps and analysis."""
import argparse,json
from fpm import dynamics,parameters,averages,heatmaps,published

def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare-data');p.add_argument('--archive',required=True);p.add_argument('--update-archive');p.add_argument('--check-only',action='store_true')
    p=sub.add_parser('build-engine');p.add_argument('--compiler',required=True)
    p=sub.add_parser('validate');p.add_argument('--design',required=True)
    p=sub.add_parser('plan');p.add_argument('--template',required=True);p.add_argument('--grid',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('simulate');p.add_argument('--design',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('averages');p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--groups',nargs='+',required=True);p.add_argument('--values',nargs='+',required=True)
    p=sub.add_parser('heatmap');p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--x',required=True);p.add_argument('--y',required=True);p.add_argument('--value',required=True)
    p=sub.add_parser('render');p.add_argument('--group',choices=['all','response','topology','sensitivity'],default='all')
    p=sub.add_parser('reproduce-published');p.add_argument('--checks-only',action='store_true')
    a=parser.parse_args()
    if a.command=='prepare-data':
        from fpm import data_bundle
        result=data_bundle.prepare(a.archive,a.check_only,a.update_archive)
    elif a.command=='build-engine':result=dynamics.build_engine(a.compiler)
    elif a.command=='validate':result=dict(rows=len(dynamics.read_design(a.design)),status='VALID')
    elif a.command=='plan':result=parameters.make_sweep(a.template,a.grid,a.output)
    elif a.command=='simulate':result=dynamics.run_design(a.design,a.output)
    elif a.command=='averages':result=averages.summarize(a.input,a.output,a.groups,a.values)
    elif a.command=='heatmap':result=heatmaps.measured_matrix(a.input,a.output,a.x,a.y,a.value)
    elif a.command=='render':
        from fpm.visuals import render
        from pathlib import Path
        receipt=render(Path(__file__).resolve().parent,a.group)
        result=dict(status=receipt['status'],group=a.group,figure_assets=receipt['figure_assets'])
    else:return published.reproduce(a.checks_only)
    print(json.dumps(result));return 0

if __name__=='__main__':raise SystemExit(main())
