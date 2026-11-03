"""Directed numerical checks on retained results; no epidemic simulation."""
from pathlib import Path
import hashlib,json,re
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1]
def save(name,x):
    p=R/name;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def table(name,x):
    p=R/name;p.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(x).to_csv(p,index=False,float_format='%.17g')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
raw=pd.read_csv(R/'evidence/fig4/source/M27_G2_RUN_LEVEL_DERIVED.tsv',sep='\t',low_memory=False)
assert len(raw)==45056 and raw.run_fingerprint.nunique()==45056
assert raw.censored.eq(0).all() and raw.stop_reason.eq('absorbed_I0').all()
assert np.array_equal(raw.ever_infected_count,raw.initial_i+raw.infections_into_n+raw.infections_into_q)
assert np.array_equal(raw.outbreak_size,raw.ever_infected_count)
for derived,counts in [('reach','ever_infected_count'),('imax_fraction','peak_infected'),('q_final_fraction','final_q')]:
    assert np.allclose(raw[derived],raw[counts]/raw.population,rtol=0,atol=1e-12)
assert np.allclose(raw.q_time_average,raw.q_person_time/raw.duration,rtol=0,atol=1e-12)
assert ((raw.q_person_time>=0)&(raw.q_person_time<=raw.duration+1e-8)).all()
V=pd.read_csv(R/'reproduction/main/data/fig4_values.csv');display=V[V.within_historical_domain].copy()
col={'reach':'reach','peak':'imax_fraction','duration':'duration','q_final':'q_final_fraction'}
errs=[]
for (top,beta),g in raw.groupby(['topology','beta_n']):
    v=V[(V.system==top)&np.isclose(V.beta_n,beta)]
    assert len(g)==128 and len(v)==4
    assert g.graph_sha256.nunique()==(8 if top.startswith('rr_') else 1)
    for obs,c in col.items():errs.append(abs(g[c].mean()-v[v.observable==obs]['mean'].iloc[0]))
assert max(errs)<1e-10
wm=pd.read_csv(R/'evidence/fig4/source/M28_E6_WELL_MIXED_P5_COORDS_SUMMARY.tsv',sep='\t',low_memory=False)
wm=wm[np.isclose(wm.delta,10)&np.isclose(wm.q0,.5)]
assert len(wm)==720 and wm.censored.eq(0).all()
werr=[]
for beta,g in wm.groupby('beta_n'):
    assert len(g)==16
    for obs,c in [('reach','ever_infected_count'),('peak','peak_infected'),('duration','duration'),('q_final','final_q')]:
        x=g[c].to_numpy(float)/(g.population.to_numpy(float) if obs!='duration' else 1)
        t=V[(V.system=='well_mixed')&(V.observable==obs)&np.isclose(V.beta_n,beta)].iloc[0]
        werr.extend([abs(x.mean()-t['mean']),abs(x.mean()-1.96*x.std(ddof=1)/4-t.ci_low),abs(x.mean()+1.96*x.std(ddof=1)/4-t.ci_high)])
assert max(werr)<1e-10
bd=float((R/'reproduction/main/data/fig4_beta_dagger.txt').read_text())
display['drawn_in_article']=True
display.loc[(display.system=='ode')&(display.observable!='peak')&(display.beta_n>bd-.05)&(display.beta_n<bd+.05),'drawn_in_article']=False
table('inventory/FIG4_EXACT_DRAWING_DOMAIN.csv',display)
slopes=[]
for (system,obs),g in display.groupby(['system','observable']):
    g=g.sort_values('beta_n');x=g.beta_n.to_numpy();y=g['mean'].to_numpy()
    for j in range(1,len(g)):
        slopes.append(dict(system=system,observable=obs,beta_left=x[j-1],beta_right=x[j],delta_beta=x[j]-x[j-1],left=y[j-1],right=y[j],delta_mean=y[j]-y[j-1],slope=(y[j]-y[j-1])/(x[j]-x[j-1]),both_drawn=bool(g.drawn_in_article.iloc[j-1] and g.drawn_in_article.iloc[j])))
table('diagnostics/FIG4_ALL_ADJACENT_VALUES.csv',slopes)
s=pd.DataFrame(slopes)
top=s.assign(abs_slope=s.slope.abs()).sort_values('abs_slope',ascending=False).groupby(['system','observable']).head(1)
table('diagnostics/FIG4_LARGEST_MEASURED_SLOPE_PER_SERIES.csv',top)
ode=V[V.system=='ode'].pivot(index='beta_n',columns='observable',values='mean')
jump=[dict(observable=o,beta_left=4.4,beta_right=4.45,left=float(ode.loc[4.4,o]),right=float(ode.loc[4.45,o]),difference=float(ode.loc[4.45,o]-ode.loc[4.4,o])) for o in col]
table('diagnostics/FIG4_ODE_BRANCH_CHANGE.csv',jump)
grids=[tuple(g.beta_n.sort_values()) for _,g in display[~display.system.isin(['ode','well_mixed'])].groupby(['system','observable'])]
assert len(grids)==32 and len(set(grids))==1 and len(grids[0])==39
assert np.isfinite(display['mean']).all() and not display.duplicated(['system','observable','beta_n']).any()
finite=display[display.system!='ode'];assert finite[['ci_low','ci_high']].notna().all().all()
assert ((finite.ci_low<=finite['mean'])&(finite['mean']<=finite.ci_high)).all()
save('qa/FIG4_DIRECTED_CHECK.json',dict(status='PASS_WITH_ENDPOINT_AND_DESIGN_LIMITS',network_runs=len(raw),display_network_runs=int(raw.beta_n.le(20).sum()),network_cells=352,display_network_cells=312,network_mean_checks=len(errs),max_network_mean_error=max(errs),wm_runs=720,wm_cells=45,max_wm_mean_and_CI_error=max(werr),ODE_source_coordinates=381,drawn_records=len(display[display.drawn_in_article]),withheld_ODE_records=len(display[~display.drawn_in_article]),withheld_beta=[4.4,4.45],affected_panels=['A','C','D'],branch_reference=bd,network_grid_complete=True,network_grid_points_per_series=39,ODE_peak_continuous_drawing=True,new_epidemic_runs=0,new_ODE_integrations=0,prior_ODE_reintegration='9 anchors, 36 observables; copied evidence only, not rerun in this task',interpretation='Jump is a change in the first-passage observation functional, not a discontinuous ODE state trajectory. Adjacent finite means are measured coordinates; joining lines are guides. Different N and endpoints do not isolate a pure topology effect.'))
# Point estimates in Fig11, including all five contrasts, use the same M27 ledger.
c=pd.read_csv(R/'reproduction/main/data/fig11_contrasts.csv');cerr=[]
means=raw.groupby(['topology','beta_n']).reach.mean()
for t in c.itertuples():cerr.append(abs(means.loc[(t.spatial,t.beta_n)]-means.loc[(t.rrg,t.beta_n)]-t.difference))
assert len(c)==220 and max(cerr)<1e-10
save('qa/M27_SHARED_OBSERVABLES.json',dict(status='PASS',unique_runs=len(raw),reach_count_identity_checks=len(raw),protection_person_time_bounds_checks=len(raw),derived_variable_checks=4*len(raw),Fig11_contrasts_checked=220,Fig11_contrasts_drawn=195,max_contrast_error=max(cerr),calibrated_new_graph_bootstrap=False))
# Source-level checks document the actual simulator; they are not binary attestation.
src=R/'engine/src/fpm_network_sim.cpp';txt=src.read_text()
checks=[('stop_at_I0','while (infected.size() > 0'),('physical_time_protection_integral','double(count[SQ] + count[IQ]) / graph.n * interval'),('single_shuffle','signal_rng.shuffle(signal_donor)'),('no_focal_in_BFS','seen[root] = root;'),('infection_count','ever_infected_count++')]
entries=[]
for name,needle in checks:
    assert needle in txt;line=txt[:txt.index(needle)].count('\n')+1
    entries.append(dict(check=name,path=src.relative_to(R).as_posix(),line=line,matched=needle,scope='source inspection; not a new simulator validation'))
save('qa/ENGINE_DIRECTED_SOURCE_CHECKS.json',entries)

print("Shared-observable checks PASS; no epidemic runs.")
