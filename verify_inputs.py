"""Verify selected scientific inputs without opening historical source paths."""
from pathlib import Path
import argparse,csv,hashlib,json
R=Path(__file__).resolve().parent

def main(include_generated=False):
    rows=list(csv.DictReader((R/'provenance/SELECTED_INPUTS.csv').open(encoding='utf-8')))
    failed=[];checked=0;deferred=0
    for row in rows:
        generated=row['path'].startswith(('reproduction/main/manuscript/figures/','reproduction/si/supplement/figures/'))
        if generated and not include_generated:
            deferred+=1;continue
        p=R/row['path']
        checked+=1
        if not p.is_file() or p.stat().st_size!=int(row['bytes']) or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:failed.append(row['path'])
    assert not failed,failed
    print(json.dumps(dict(status='PASS_SELECTED_INPUTS',files=checked,generated_outputs_deferred=deferred,historical_path_reads=False)))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--include-generated',action='store_true')
    main(parser.parse_args().include_generated)
