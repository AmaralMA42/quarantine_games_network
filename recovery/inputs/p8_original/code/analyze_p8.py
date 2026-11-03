#!/usr/bin/env python3
"""Read-only scientific analysis of verified P8 RRG3 output; graph -> replicate paired hierarchical bootstrap.
No CTMC invocation. Intended to execute reproducibly with original 16 shards.
"""
import csv,glob,hashlib,json,math,pathlib,sys
import numpy as np
import pandas as pd
BASE=pathlib.Path(__file__).resolve().parents[1];E=BASE/'execution';A=BASE/'analysis';A.mkdir(exist_ok=True)
S=BASE/'design/P8_256_FULL_DESIGN.tsv'
rs=pd.concat([pd.read_csv(p,sep='\t',dtype={'run_id':str,'run_fingerprint':str,'graph_sha256':str}) for p in sorted(E.glob('P8*_result.tsv'))],ignore_index=True)
ds=pd.read_csv(S,sep='\t',dtype=str)
assert len(rs)==len(ds)==256 and rs['censored'].astype(int).eq(0).all() and rs.stop_reason.eq('absorbed_I0').all()
keys=['run_fingerprint','run_id','population','graph_sha256','q_init_seed','i_init_seed','dynamic_seed','signal_seed','beta_n','behavior_mode']
r=rs.merge(ds[keys].rename(columns={c:c+'_design' for c in keys if c not in ('run_id','run_fingerprint')}),on=['run_id','run_fingerprint'],validate='one_to_one')
assert len(r)==256
for k in keys[2:]:assert all(r[k].astype(str)==r[k+'_design'].astype(str)),k
assert r['initial_i'].astype(int).eq(1).all()
assert r['initial_q'].astype(int).eq(16384).all()
for k in ['q_to_n_events','n_to_q_events']:assert r.loc[r.behavior_mode=='frozen',k].astype(int).eq(0).all()
r['R']=r.ever_infected_count/r.population
r['peak_frac']=r.peak_infected/r.population
r['graph_no']=r.run_id.str.extract(r'_g(\d\d)_')[0].astype(int)
r['rep_no']=r.run_id.str.extract(r'_r(\d\d)_')[0].astype(int)
assert not r[['graph_no','rep_no','beta_n','behavior_mode']].duplicated().any()
assert set(r.beta_n)=={20,50};assert set(r.behavior_mode)=={'adaptive','frozen'}
seed_keys=['q_init_seed','i_init_seed','dynamic_seed','signal_seed']
for _,gr in r.groupby(['graph_no','rep_no']):
 assert len(gr)==4
 for k in seed_keys:assert gr[k].nunique()==1
allrows=r.sort_values(['graph_no','rep_no','beta_n','behavior_mode']).copy()
allrows.to_csv(A/'P8_256_RUNLEVEL_AUDITED.csv',index=False,float_format='%.17g',lineterminator='\n')
wide=r.pivot(index=['graph_no','rep_no'],columns=['beta_n','behavior_mode'],values='R')
assert len(wide)==64
# 4 by 16: 0=B20 Adaptive,1=B20 Frozen,2=B50 Adaptive,3=B50 Frozen
m=np.zeros((4,16,4));pairs=[]
for g in range(4):
 for j in range(16):
  w=wide.loc[(g,j)]
  aa=float(w[(20,'adaptive')]);af=float(w[(20,'frozen')]);ba=float(w[(50,'adaptive')]);bf=float(w[(50,'frozen')]);m[g,j]=[aa,af,ba,bf]
  pairs.append({'graph':g,'replicate':j,'beta20_adaptive':aa,'beta20_frozen':af,'D20_A_minus_F':aa-af,'beta50_adaptive':ba,'beta50_frozen':bf,'D50_A_minus_F':ba-bf,'interaction_D50_minus_D20':(ba-bf)-(aa-af)})
pd.DataFrame(pairs).to_csv(A/'P8_64_CROSSED_SEED_PAIRS.csv',index=False,float_format='%.17g',lineterminator='\n')
vg=m.mean(axis=1);gtable=[]
for g in range(4):
 v=vg[g];gtable.append({'graph':g,'sha256':r.loc[r.graph_no==g,'graph_sha256'].iloc[0],'R_A20':v[0],'R_F20':v[1],'D20':v[0]-v[1],'R_A50':v[2],'R_F50':v[3],'D50':v[2]-v[3],'interaction':(v[2]-v[3])-(v[0]-v[1]),'pairs_per_beta':16})
pd.DataFrame(gtable).to_csv(A/'P8_FOUR_GRAPH_EFFECTS.csv',index=False,float_format='%.17g',lineterminator='\n')
# Graph-wise paired bootstrap hierarchical; same g and replicate index used across both betas
B=10000; rng=np.random.default_rng(20261012);boots=np.empty((B,7))
for ib in range(B):
 gs=rng.integers(0,4,4)
 samp=np.stack([m[g,rng.integers(0,16,16)] for g in gs])
 v=samp.mean(axis=(0,1));boots[ib]=[v[0],v[1],v[2],v[3],v[0]-v[1],v[2]-v[3],(v[2]-v[3])-(v[0]-v[1])]
terms=['R_A20','R_F20','R_A50','R_F50','D20','D50','interaction'];observed=m.mean(axis=(0,1));centers=list(observed)+[float(observed[0]-observed[1]),float(observed[2]-observed[3]),float((observed[2]-observed[3])-(observed[0]-observed[1]))];st=[]
for i,t in enumerate(terms):
 lo,hi=np.quantile(boots[:,i],[.025,.975]);st.append({'estimand':t,'point':centers[i],'ci95_lower_graph_boot':float(lo),'ci95_upper_graph_boot':float(hi),'n_graphs':4,'n_pairs_per_graph_beta':16,'bootstrap_B':B,'fixed_seed':20261012,'CI_descriptive_exploratory':True})
pd.DataFrame(st).to_csv(A/'P8_HIERARCHICAL_BOOTSTRAP_SUMMARY.csv',index=False,float_format='%.17g',lineterminator='\n')
conditions=r.groupby(['beta_n','behavior_mode']).agg(n=('R','size'),mean_reach=('R','mean'),median_reach=('R','median'),reach_02=('R',lambda x:int(sum(x>=.02))),reach_05=('R',lambda x:int(sum(x>=.05))),reach_10=('R',lambda x:int(sum(x>=.1))),duration_mean=('duration','mean'),peak_mean=('peak_frac','mean'),final_q_mean=('final_q','mean'),max_q_mean=('max_q','mean'),q_to_n_mean=('q_to_n_events','mean'),n_to_q_mean=('n_to_q_events','mean')).reset_index()
conditions.to_csv(A/'P8_ARMS_AND_DETECTORS.csv',index=False,float_format='%.17g',lineterminator='\n')
qa={'state':'P8_256_COMPLETE_VERIFIED_UNPROMOTED','n_runs':len(r),'n_unique_ids':r.run_id.nunique(),'n_unique_fp':r.run_fingerprint.nunique(),'n_graphs':4,'n_pairs_g_beta':16,'n_pairs_beta':64,'n_censored':int(r.censored.astype(int).sum()),'n_absorbed':int(sum(r.stop_reason=='absorbed_I0')),'all_paired_four_seeds':True,'source_graph_archive_sha256':'d3302527263d403e2a0de9b0685b50eedcc3fbf3b8e00717d5f792377b951a66','binary_sha256':'3fe1a6c82910632597ac8574ab94ea0e94eb9aeceedd3bda39d5afa837b8194e','graph_seed_original_available':False,'graph_seed_field':0,'bounds':{'D20':centers[4],'D50':centers[5],'interaction':centers[6]},'source_cross_checks':'design/output 256/256 and FROZEN Q/N 0/128','bootstrap':'graph->pair, 10000 draws, 4 graph clusters; across beta uses same graph and replicate indices','rc94_unchanged':True,'g3_confirmatory_not_started':True}
(A/'P8_INDEPENDENT_QA.json').write_text(json.dumps(qa,indent=2,sort_keys=True)+'\n')
print('P8',json.dumps(qa,indent=2));print(pd.DataFrame(st).to_string(index=False));print(pd.DataFrame(gtable).to_string(index=False))
