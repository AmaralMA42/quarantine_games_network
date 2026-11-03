"""Measured-data reconstruction through shared modules and numerical audits."""
from pathlib import Path
import json,subprocess,sys

def reproduce(checks_only=False):
    root=Path(__file__).resolve().parents[1]
    for path in ['qa','inventory','diagnostics','tmp','recovery/qa','recovery/results','reproduction/main/qa','reproduction/main/manuscript/figures','reproduction/si/qa','reproduction/si/supplement/figures']:(root/path).mkdir(parents=True,exist_ok=True)
    def check(script,args=()):
        result=subprocess.run([sys.executable,str(root/script)]+list(args),cwd=root,capture_output=True,text=True,errors='replace',timeout=900)
        (root/'qa'/(Path(script).stem+'_STDOUT.txt')).write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
        assert result.returncode==0,(script,result.stdout[-2000:],result.stderr[-3000:])
    check('verify_inputs.py')
    if not checks_only:
        from .visuals import render
        figures=render(root)
    for script in ['analysis/check_observables.py','analysis/check_recovered.py']:check(script)
    check('verify_inputs.py',['--include-generated'] if not checks_only else [])
    receipt=dict(status='PASS_UNIFIED_DATA_REPLAY',figure_assets=figures['figure_assets'] if not checks_only else 'not requested',scientific_inputs_unchanged=True,legacy_renderer_execution=False,legacy_renderer_import=False,new_article_epidemic_runs=0)
    (root/'qa/REPRODUCTION.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt));return 0
