"""Verify immutable historical observations and the current figure selection."""
from pathlib import Path
import argparse,csv,hashlib,json
R=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(include_generated=False):
    rows=list(csv.DictReader((R/'provenance/SELECTED_INPUTS.csv').open(encoding='utf-8')))
    failed=[];checked=0;historical_generated=0
    for row in rows:
        generated=row['path'].startswith(('reproduction/main/manuscript/figures/','reproduction/si/supplement/figures/'))
        if generated:
            historical_generated+=1;continue
        p=R/row['path'];checked+=1
        if not p.is_file() or p.stat().st_size!=int(row['bytes']) or sha(p)!=row['sha256']:failed.append(row['path'])
    config=json.loads((R/'ARTICLE_UPDATE.json').read_text('utf-8'))
    manifest=R/config['manifest_path']
    assert manifest.is_file(),'Attach both the base data and the current figure update.'
    assert sha(manifest)==config['manifest_sha256']
    current=list(csv.DictReader(manifest.open(encoding='utf-8',newline='')))
    for row in current:
        p=R/row['path']
        assert p.resolve().is_relative_to(R.resolve()),row['path']
        if not p.is_file() or p.stat().st_size!=int(row['bytes']) or sha(p)!=row['sha256']:failed.append(row['path'])
    active=json.loads((R/'article_update/ACTIVE_FIGURE_ASSETS.json').read_text('utf-8'))
    if include_generated:
        for row in active:
            p=R/row['generated']
            if not p.is_file() or p.stat().st_size!=row['bytes'] or sha(p)!=row['sha256']:failed.append(row['generated'])
    assert not failed,failed
    print(json.dumps(dict(status='PASS_SELECTED_INPUTS',historical_scientific_inputs=checked,
        historical_generated_outputs_replaced_in_current_selection=historical_generated,
        current_input_and_reference_files=len(current),current_generated_outputs_checked=len(active) if include_generated else 0,
        original_provenance_manifest_rewritten=False,historical_path_reads=False)))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--include-generated',action='store_true')
    main(parser.parse_args().include_generated)
