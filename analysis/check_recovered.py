"""Reanalyze retained individual records and original graph identities."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
S=R/'recovery/inputs/p8_original';W=R/'tmp/p8_reanalysis';W.mkdir(parents=True,exist_ok=True)
for folder in ['code','execution','design']:shutil.copytree(S/folder,W/folder,dirs_exist_ok=True)
proc=subprocess.run([sys.executable,str(W/'code/analyze_p8.py')],capture_output=True,text=True,errors='replace',timeout=180)
assert proc.returncode==0,proc.stderr[-2000:]
products=[]
for source in sorted((S/'analysis').glob('*')):
    if source.suffix not in {'.csv','.json'}:continue
    target=W/'analysis'/source.name
    same=source.read_bytes()==target.read_bytes()
    semantic=json.loads(source.read_text())==json.loads(target.read_text()) if source.suffix=='.json' else same
    assert semantic,source.name
    products.append(dict(file=source.name,bytes_equal=same,semantic_equal=semantic))
for name in ['P8_FOUR_GRAPH_EFFECTS.csv','P8_HIERARCHICAL_BOOTSTRAP_SUMMARY.csv']:
    assert sha(S/'analysis'/name)==sha(R/'reproduction/si/data/single_introduction'/name)
raw=pd.read_csv(R/'recovery/inputs/s4_original/analysis/I1_M26_PLUS_NEW_832_RUN_LEVEL.csv',dtype=str)
assert len(raw)==raw.run_fingerprint.nunique()==832
for c in ['beta_n','ever_infected_count','population','initial_i','censored']:raw[c]=pd.to_numeric(raw[c])
assert raw.initial_i.eq(1).all() and raw.population.eq(33120).all()
assert raw.censored.eq(0).all() and raw.stop_reason.eq('absorbed_I0').all()
raw['reach']=raw.ever_infected_count/raw.population
reference=pd.read_csv(R/'reproduction/si/data/s4/PRIMARY_832_COUNTS.csv');checks=0
for row in reference.itertuples():
    d=raw[np.isclose(raw.beta_n,row.beta_n)];assert len(d)==row.n
    for cutoff,key in [(.02,'n_ge2pct'),(.05,'n_ge5pct'),(.10,'n_ge10pct')]:
        assert int(d.reach.ge(cutoff).sum())==getattr(row,key);checks+=1
proc=subprocess.run([sys.executable,str(R/'recovery/regenerate_fig10_graphs.py')],capture_output=True,text=True,errors='replace',timeout=180)
assert proc.returncode==0,proc.stderr[-2000:]
graphs=json.loads((R/'recovery/qa/FIG10_GENERATOR_REPLAY.json').read_text())
for item in graphs['files']:
    assert sha(R/'recovery/results/fig10_graphs'/item['file'])==item['sha256']
report=dict(status='PASS_RECOVERED_RECORDS',P8_individual_outcomes=256,P8_matched_pairs=128,P8_graphs=4,P8_products=products,P8_graph_generator_seeds='not recorded',S4_individual_outcomes=832,S4_cutoff_count_checks=checks,S4_dense_time_traces_claimed=False,Fig10_original_graphs_regenerated=12,new_epidemic_runs=0)
(R/'qa/RECOVERED_RECORDS.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='P8_products'}))
