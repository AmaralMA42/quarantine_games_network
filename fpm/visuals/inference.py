from .style import ORDER as NETWORK_ORDER, STYLE as NETWORK_STYLE
from pathlib import Path
import json, hashlib, math
import numpy as np, pandas as pd
from scipy.stats import beta, binom, binomtest
from scipy.optimize import brentq

def conditional_establishment(root):
    R = Path(root) / 'reproduction/main/preserved_renderers/fig3'
    ORDER = NETWORK_ORDER
    STYLE = NETWORK_STYLE
    Q = R / 'qa'
    Q.mkdir(exist_ok=True)
    S = R / 'results'
    S.mkdir(exist_ok=True)
    for b in json.loads((R / 'provenance/INPUT_BINDINGS.json').read_text()):
        assert hashlib.sha256((R / b['stored_path']).read_bytes()).hexdigest() == b['sha256']
    cols = ['run_id', 'run_fingerprint', 'topology', 'beta_n', 'graph_sha256', 'population', 'ever_infected_count', 'final_q', 'censored', 'stop_reason']
    m = pd.read_csv(R / 'inputs/m27/runs.tsv', sep='\t', usecols=cols)
    m['campaign'] = 'M27'
    parts = []
    for p in sorted((R / 'inputs/ab').glob('*.tsv')):
        d = pd.read_csv(p, sep='\t', usecols=cols)
        d['campaign'] = 'AB'
        parts.append(d)
    a = pd.concat(parts, ignore_index=True)
    D = pd.concat([m, a], ignore_index=True)
    assert len(m) == 45056 and len(a) == 2560 and (D.run_fingerprint.nunique() == 47616) and (D.run_id.nunique() == 47616)
    assert D.censored.eq(0).all() and D.stop_reason.eq('absorbed_I0').all()
    D['cell_id'] = [f'{c}:{t}:{b:.8f}' for c, t, b in zip(D.campaign, D.topology, D.beta_n)]
    D['established'] = (20 * D.ever_infected_count >= D.population).astype(int)
    D['reach'] = D.ever_infected_count / D.population
    D['q_final'] = D.final_q / D.population
    E = pd.read_csv(R / 'data/FIG3_P_MACRO_392_MEASURED_POINTS.csv')
    H = pd.read_csv(R / 'data/FIG3_REACH_392_MEASURED_POINTS.csv').set_index('cell_id')
    F = pd.read_csv(R / 'data/FIG3_Q_FINAL_392_MEASURED_POINTS.csv').set_index('cell_id')
    assert len(E) == 392 and set(D.cell_id) == set(E.cell_id)

    def cp(x, n):
        return (0.0 if x == 0 else float(beta.ppf(0.025, x, n - x + 1)), 1.0 if x == n else float(beta.ppf(0.975, x + 1, n - x)))

    def kl(q, n):
        c = math.log(40) / n
        if q == 0:
            return (0.0, -math.expm1(-c))
        if q == 1:
            return (math.exp(-c), 1.0)

        def fun(p):
            return q * math.log(q / p) + (1 - q) * math.log((1 - q) / (1 - p)) - c
        return (brentq(fun, 1e-300, q, xtol=1e-14), brentq(fun, q, 1 - np.finfo(float).eps, xtol=1e-14) if fun(1 - np.finfo(float).eps) > 0 else 1.0)
    graphs = []
    out = []
    maxmean = 0.0
    reg = []
    for e in E.itertuples():
        d = D[D.cell_id == e.cell_id]
        groups = list(d.groupby('graph_sha256', sort=True))
        ng = len(groups)
        n = len(d)
        x = int(d.established.sum())
        q = x / n
        assert ng == e.n_graphs and n == e.n and (len(set((len(g) for _, g in groups))) == 1)
        assert all((len(g) == e.dynamics_per_graph for _, g in groups))
        assert ng == (8 if e.topology.startswith('rr_') else 1)
        assert abs(q - e.estimate) < 1e-14
        er = abs(d.reach.mean() - H.loc[e.cell_id, 'estimate'])
        ef = abs(d.q_final.mean() - F.loc[e.cell_id, 'estimate'])
        maxmean = max(maxmean, er, ef)
        assert max(er, ef) < 1e-14
        for h, g in groups:
            graphs.append(dict(cell_id=e.cell_id, campaign=e.campaign, topology=e.topology, beta_n=e.beta_n, graph_sha256=h, n=len(g), successes=int(g.established.sum()), probability=float(g.established.mean())))
        bounds = cp(x, n) if ng == 1 else kl(q, n)
        if ng == 1:
            check = binomtest(x, n).proportion_ci(confidence_level=0.95, method='exact')
            assert max(abs(check.low - bounds[0]), abs(check.high - bounds[1])) < 1e-09
        assert 0 <= bounds[0] <= q + 1e-14 and q <= bounds[1] + 1e-14 <= 1 + 1e-14 and (bounds[1] > bounds[0])
        row = e._asdict()
        row.update(successes=x, conditional_estimate=q, conditional_lo=bounds[0], conditional_hi=bounds[1], conditional_method='EXACT_CLOPPER_PEARSON_FIXED_SPATIAL_GRAPH' if ng == 1 else 'HETEROGENEOUS_BERNOULLI_CHERNOFF_KL_FIXED_RRG_BANK', conditional_target='Equal-weight trajectory probability over the observed graph bank' if ng > 1 else 'Trajectory probability of the observed spatial graph', pointwise_level=0.95, graph_population_interval=False)
        out.append(row)
        reg.append(dict(cell_id=e.cell_id, p_macro_error=abs(q - e.estimate), reach_error=er, q_final_error=ef))

    def save(name, x):
        pd.DataFrame(x).to_csv(S / name, index=False, float_format='%.17g')
    P = pd.DataFrame(out)
    P.to_csv(R / 'data/FIG3_ESTABLISHMENT_CONDITIONAL_392.csv', index=False, float_format='%.17g')
    save('GRAPH_LEVEL_ESTABLISHMENT_COUNTS.csv', graphs)
    save('ALL_392_POINT_ESTIMATE_REGRESSION.csv', reg)
    D.to_csv(S / 'RUN_47616_OUTCOME_BINDINGS.csv', index=False, float_format='%.17g')
    cases = []
    vectors = [np.zeros(8), np.ones(8)] + [np.full(8, p) for p in [0.001, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99]] + [np.linspace(0, 1, 8), np.array([0, 0, 0, 0, 1, 1, 1, 1]), np.array([0.001, 0.001, 0.01, 0.01, 0.1, 0.1, 0.5, 0.99])]
    for per in [8, 16]:
        n = 8 * per
        intervals = np.array([kl(k / n, n) for k in range(n + 1)])
        for index, ps in enumerate(vectors):
            dist = np.array([1.0])
            for p in ps:
                dist = np.convolve(dist, binom.pmf(np.arange(per + 1), per, p))
            mu = float(np.mean(ps))
            coverage = float(dist[(intervals[:, 0] <= mu + 1e-12) & (intervals[:, 1] >= mu - 1e-12)].sum())
            assert abs(dist.sum() - 1) < 1e-10 and coverage >= 0.95 - 1e-10
            cases.append(dict(method='RRG_CONDITIONAL_KL', graphs=8, per_graph_n=per, case=index, probabilities=';'.join((f'{p:g}' for p in ps)), true_fixed_bank_mean=mu, exact_coverage=coverage, label='HYPOTHETICAL_PROBABILITY_SUM_NOT_OBSERVED_FPM_DATA'))
    for n in [64, 128]:
        ints = np.array([cp(k, n) for k in range(n + 1)])
        minimum = 1.0
        for p in np.linspace(1e-05, 0.99999, 1001):
            coverage = binom.pmf(np.arange(n + 1), n, p)[(ints[:, 0] <= p) & (ints[:, 1] >= p)].sum()
            minimum = min(minimum, float(coverage))
            assert coverage >= 0.95 - 1e-10
        cases.append(dict(method='SPATIAL_EXACT_CP', per_graph_n=n, case='1001_probability_grid', exact_coverage=minimum, label='MINIMUM_ON_HYPOTHETICAL_GRID_NOT_FPM_COVERAGE_ESTIMATE'))
    save('HYPOTHETICAL_EXACT_COVERAGE_DIAGNOSTICS.csv', cases)
    visible = P[P.beta_n <= 2.6]
    assert len(visible) == 192
    report = dict(status='PASS_CONDITIONAL_METHOD', runs=47616, cells=392, graph_rows=len(graphs), spatial_cells=int((P.n_graphs == 1).sum()), rrg_cells=int((P.n_graphs == 8).sum()), displayed_establishment_points=192, displayed_spatial_points=120, displayed_rrg_points=72, displayed_C_D_mean_points=672, central_estimates_unchanged=True, max_continuous_mean_error=maxmean, zero_width_selected_establishment_intervals=0, conditional_on_observed_graphs=True, independent_conditional_trajectories_assumed=True, coverage_over_new_graphs_claimed=False, coverage_simultaneous_claimed=False, hypothetical_exact_diagnostic_cases=len(cases), new_epidemic_runs=0, hypothetical_observations_added_to_data=0, C_D_uncertainty_reestimated=False, source_C_D_intervals_retained=True)
    (Q / 'CONDITIONAL_INFERENCE_QA.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
