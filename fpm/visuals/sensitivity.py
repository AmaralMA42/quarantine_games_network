"""Scientific compositions using shared identities and export primitives."""
from pathlib import Path
import os, json, hashlib
import numpy as np, pandas as pd, matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.container import ErrorbarContainer
from matplotlib.patches import Rectangle, FancyBboxPatch
from PIL import Image
from scipy.signal import find_peaks
from statistics import NormalDist
from .style import STATE_COLORS, LANDMARKS

def endpoint_sensitivity(ctx):
    R = ctx.root
    ROOT = ctx.root
    O = ctx.output
    P = ctx.preview
    ORDER = ctx.order
    STYLE = ctx.style
    NAMES = [STYLE[t][0] for t in ORDER]
    qa = ctx.qa
    coordinates = ctx.coordinates
    clean = ctx.clean
    handles = ctx.handles
    legend = ctx.legend
    panel = ctx.panel
    save = ctx.save
    topology_legend = ctx.topology_legend
    d = pd.read_csv(R / 'data/s1/ENDPOINT_THRESHOLD_VALUES.csv')
    assert len(d) == 123 and d.N.eq(16384).all() and d['hit'].all()
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.4))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.23, top=0.88, wspace=0.37)
    records = []
    for factor, color, ls in [(0.5, '#999999', ':'), (1, '#222222', '-'), (2, '#666666', '--')]:
        t = d[d.factor == factor].sort_values('beta_n')
        assert len(t) == 41
        for ax, obs, letter, label in zip(axs, ['T_epi', 'reach', 'q'], 'ABC', ['$T_{\\rm epi}$', '$R_\\dagger$', '$q_\\dagger$']):
            ax.plot(t.beta_n, t[obs], marker='o', ms=2, markeredgewidth=0, linestyle=ls, color=color, lw=0.8, label=f'$c={factor:g}$')
            ax.set(xlim=(3.5, 5.5), xticks=[3.5, 4, 4.5, 5, 5.5], xlabel='$\\beta_N$', ylabel=label)
            records.extend((dict(panel=letter, observable=obs, beta_n=float(x), factor=factor, value=float(y)) for x, y in zip(t.beta_n, t[obs])))
    leg = axs[0].legend(loc='center right', bbox_to_anchor=(0.99, 0.53), frameon=True, fancybox=False, framealpha=1, edgecolor='.7', fontsize=7, handlelength=1, handletextpad=0.3, borderpad=0.3, labelspacing=0.35)
    fig.canvas.draw()
    box = leg.get_window_extent(fig.canvas.get_renderer())
    assert not any((line.get_path().transformed(line.get_transform()).intersects_bbox(box, filled=False) for line in axs[0].lines)), 'S1 legend overlaps a measured curve or connecting guide.'
    for ax, letter in zip(axs, 'ABC'):
        clean(ax)
        panel(ax, letter)
    save(fig, 'S1_ENDPOINT_THRESHOLD_STANDARD', records)

def mean_peak_map(ctx):
    R = ctx.root
    ROOT = ctx.root
    O = ctx.output
    P = ctx.preview
    ORDER = ctx.order
    STYLE = ctx.style
    NAMES = [STYLE[t][0] for t in ORDER]
    qa = ctx.qa
    coordinates = ctx.coordinates
    clean = ctx.clean
    handles = ctx.handles
    legend = ctx.legend
    panel = ctx.panel
    save = ctx.save
    topology_legend = ctx.topology_legend
    source = pd.read_csv(R / 'data/s2/PAPER_WINDOW_5OBS.csv')
    d = source[source.observable == 'robust_peak_count'].copy()
    assert len(d) == 4480
    fig, axs = plt.subplots(2, 4, figsize=(7.1, 3.6))
    fig.subplots_adjust(left=0.07, right=0.885, bottom=0.14, top=0.91, wspace=0.28, hspace=0.4)
    records = []
    for j, (ax, topo) in enumerate(zip(axs.flat, ORDER)):
        t = d[d.topology == topo]
        xs = sorted(t.beta_n.unique())
        ys = sorted(t.delta.unique())
        z = t.pivot(index='delta', columns='beta_n', values='mean').reindex(index=ys, columns=xs)
        assert z.shape == (20, 28) and z.notna().all().all()
        xe = np.r_[np.array(xs) - 0.25, xs[-1] + 0.25]
        ye = np.r_[np.array(ys) - 0.5, ys[-1] + 0.5]
        im = ax.imshow(z.to_numpy(), origin='lower', aspect='auto', extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap='viridis', vmin=0, vmax=4, interpolation='nearest', resample=False)
        assert im.get_interpolation() == 'nearest' and (not im.get_resample())
        ax.set(xticks=[1.5, 5, 10, 15], yticks=[1, 10, 20])
        panel(ax, chr(65 + j))
        if j >= 4:
            ax.set_xlabel('$\\beta_N$')
        if j % 4 == 0:
            ax.set_ylabel('$\\delta$')
        records.extend((dict(panel=chr(65 + j), topology=topo, beta_n=float(row.beta_n), delta=float(row.delta), mean=float(row['mean'])) for _, row in t.iterrows()))
    fig.colorbar(im, cax=fig.add_axes([0.915, 0.18, 0.018, 0.65]), label='$\\langle n_p\\rangle$', ticks=np.arange(0, 4.1, 0.5))
    save(fig, 'S2_MEAN_PEAK_COUNT_STANDARD', records)
    (R / 'qa/S2_GRID_RENDER.json').write_text(json.dumps(dict(interpolation='nearest', resample=False, panels=8, shape=[20, 28], measured_cells=4480, color_limits=[0, 4], cell_layer_generated_from_numeric_arrays=True, imported_image_backgrounds=0, reason='Measured-cell image avoids PDF viewer seams between adjacent vector polygons; text and axes remain vector.'), indent=2) + '\n', encoding='utf-8')

def control_sensitivity(ctx):
    R = ctx.root
    ROOT = ctx.root
    O = ctx.output
    P = ctx.preview
    ORDER = ctx.order
    STYLE = ctx.style
    NAMES = [STYLE[t][0] for t in ORDER]
    qa = ctx.qa
    coordinates = ctx.coordinates
    clean = ctx.clean
    handles = ctx.handles
    legend = ctx.legend
    panel = ctx.panel
    save = ctx.save
    topology_legend = ctx.topology_legend
    d = pd.read_csv(R / 'data/s3/ORIGINAL_PAIRED_CONTRASTS.csv')
    d = d[d.metric == 'reach']
    audit = pd.read_csv(R / 'data/s3/AUDITED_128_REACH_CONTRASTS.csv')
    keys = ['topology', 'anchor', 'contrast']
    joined = d.merge(audit, on=keys, suffixes=('', '_audit'), validate='one_to_one')
    assert len(joined) == 128
    for c in ['estimate', 'ci_low', 'ci_high']:
        assert np.allclose(joined[c], joined[c + '_audit'], atol=1e-12, rtol=0)
    panels = [('FROZEN_GLOBAL-ADAPTIVE_GLOBAL', '$\\Delta^{\\rm soc}R_\\infty$', (-0.65, 0.78)), ('ADAPTIVE_LOCAL_R4-ADAPTIVE_GLOBAL', '$\\Delta^{\\rm info}R_\\infty$', (-0.12, 0.75)), ('ADAPTIVE_LOCAL_R4-ADAPTIVE_LOCAL_SHUFFLED_R4', '$\\Delta^{\\rm align}R_\\infty$', (-0.18, 0.21))]
    anchors = ['LEFT_BETA8', 'VALLEY', 'REBOUND_FIRST', 'TAIL_BETA40']
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 3.1))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.19, top=0.91, wspace=0.35)
    records = []
    segments = []
    for ax, letter, (contrast, ylabel, ylim) in zip(axs, 'ABC', panels):
        for k, topo in enumerate(ORDER):
            t = d[(d.topology == topo) & (d.contrast == contrast)].set_index('anchor').loc[anchors]
            label, color, marker = STYLE[topo]
            x = np.arange(4) + (k - 3.5) * 0.065
            cont = ax.errorbar(x, t.estimate, yerr=[t.estimate - t.ci_low, t.ci_high - t.estimate], fmt=marker, linestyle='none', color=color, mfc=color, mec=color, markersize=3.5, capsize=1.8, elinewidth=0.65, capthick=0.65)
            for (anchor, row), seg in zip(t.iterrows(), cont.lines[2][0].get_segments()):
                assert np.allclose(seg[:, 1], [row.ci_low, row.ci_high], atol=1e-14, rtol=0)
                record = dict(panel=letter, topology=topo, anchor=anchor, beta_n=float(row.beta_n), estimate=float(row.estimate), ci_low=float(row.ci_low), ci_high=float(row.ci_high), n_pairs=int(row.n_pairs), graph_count=int(row.graph_count))
                records.append(record)
                segments.append(record)
        ax.axhline(0, color='.35', ls='--', lw=0.6, zorder=0)
        ax.set(xlim=(-0.48, 3.48), ylim=ylim, ylabel=ylabel, xticks=np.arange(4), xticklabels=['$\\beta_N=8$', 'Minimum', 'Increase', '$\\beta_N=40$'])
        clean(ax)
        panel(ax, letter)
    handles = [Line2D([], [], color=STYLE[t][1], marker=STYLE[t][2], ls='none', ms=4, label=STYLE[t][0]) for t in ORDER]
    indices = [0, 2, 4, 6, 1, 3, 5, 7]
    axs[0].legend(handles=[handles[i] for i in indices], loc='lower right', ncol=2, frameon=True, fancybox=False, framealpha=1, edgecolor='.7', fontsize=7, columnspacing=0.65, handletextpad=0.2, borderpad=0.4, labelspacing=0.3)
    assert len(records) == 96
    save(fig, 'S3_CONTROLS_STANDARD', records)
    (R / 'qa/S3_INTERVAL_SEGMENTS.json').write_text(json.dumps(segments, indent=2) + '\n', encoding='utf-8')

def single_introduction(ctx):
    R = ctx.root
    ROOT = ctx.root
    O = ctx.output
    P = ctx.preview
    ORDER = ctx.order
    STYLE = ctx.style
    NAMES = [STYLE[t][0] for t in ORDER]
    qa = ctx.qa
    coordinates = ctx.coordinates
    clean = ctx.clean
    handles = ctx.handles
    legend = ctx.legend
    panel = ctx.panel
    save = ctx.save
    topology_legend = ctx.topology_legend
    d = pd.read_csv(R / 'data/s4/PRIMARY_832_COUNTS.csv').sort_values('beta_n')
    assert len(d) == 8 and d.n.sum() == 832 and d.single_seed.eq(1).all()
    assert d.n.tolist() == [32, 32, 128, 128, 128, 128, 128, 128]
    assert d.n_ge5pct.tolist() == [0, 0, 9, 68, 93, 110, 86, 30]
    z = NormalDist().inv_cdf(0.975)
    prob = d.n_ge5pct / d.n
    den = 1 + z * z / d.n
    center = (prob + z * z / (2 * d.n)) / den
    half = z * np.sqrt(prob * (1 - prob) / d.n + z * z / (4 * d.n * d.n)) / den
    d['wilson_lo'] = center - half
    d['wilson_hi'] = center + half
    for scale in ['linear']:
        fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.5))
        fig.subplots_adjust(left=0.085, right=0.985, bottom=0.22, top=0.89, wspace=0.3)
        color = '#6F3FA0'
        records = []
        art = axs[0].errorbar(d.beta_n, prob, yerr=[prob - d.wilson_lo, d.wilson_hi - prob], fmt='8-', color=color, mfc=color, mec=color, markersize=4, capsize=2, lw=0.8, elinewidth=0.75, capthick=0.7)
        for (_, row), seg in zip(d.iterrows(), art.lines[2][0].get_segments()):
            assert np.allclose(seg[:, 1], [row.wilson_lo, row.wilson_hi], atol=1e-14, rtol=0)
            records.append(dict(panel='A', beta_n=float(row.beta_n), cutoff=0.05, n=int(row.n), successes=int(row.n_ge5pct), estimate=float(row.n_ge5pct / row.n), ci_low=float(row.wilson_lo), ci_high=float(row.wilson_hi)))
        for pct, color, ls in [(2, '#A88BC4', '--'), (5, '#6F3FA0', '-'), (10, '#43265F', ':')]:
            value = d[f'n_ge{pct}pct'] / d.n
            assert np.allclose(value, d[f'p_ge{pct}pct'], atol=5e-12, rtol=0)
            axs[1].plot(d.beta_n, value, color=color, marker='8', ms=3, linestyle=ls, lw=0.8, label=f'$c={pct / 100:g}$')
            records.extend((dict(panel='B', beta_n=float(row.beta_n), cutoff=pct / 100, n=int(row.n), successes=int(row[f'n_ge{pct}pct']), estimate=float(v), ci_low=None, ci_high=None) for (_, row), v in zip(d.iterrows(), value)))
        axs[0].set_ylabel('$P_{\\rm macro}$')
        axs[1].set_ylabel('$P(R_\\infty\\geq c)$')
        for ax, letter in zip(axs, 'AB'):
            if scale == 'log2':
                ax.set_xscale('log', base=2)
                ax.set_xticks(d.beta_n, [f'{b:g}' for b in d.beta_n])
                ax.minorticks_off()
                ax.set_xlim(1.8, 55)
            else:
                ax.set(xlim=(0, 52), xticks=[2, 8, 20, 32, 50])
            ax.set(xlabel='$\\beta_N$', ylim=(-0.035, 1.035), yticks=[0, 0.2, 0.4, 0.6, 0.8, 1])
            clean(ax)
            panel(ax, letter)
        legpos = dict(loc='upper left') if scale == 'log2' else dict(loc='lower center', bbox_to_anchor=(0.45, 0.02))
        leg = axs[1].legend(**legpos, frameon=True, fancybox=False, framealpha=1, edgecolor='.7', fontsize=7)
        fig.canvas.draw()
        box = leg.get_window_extent(fig.canvas.get_renderer())
        assert not any((line.get_path().transformed(line.get_transform()).intersects_bbox(box, filled=False) for line in axs[1].lines)), ('S4 legend overlaps a curve', scale)
        save(fig, 'S4_SINGLE_INTRO_' + scale.upper(), records, folder='supplement/figures' if scale == 'linear' else 'review/alternatives')

def single_introduction_controls(ctx):
    R = ctx.root
    ROOT = ctx.root
    O = ctx.output
    P = ctx.preview
    ORDER = ctx.order
    STYLE = ctx.style
    NAMES = [STYLE[t][0] for t in ORDER]
    qa = ctx.qa
    coordinates = ctx.coordinates
    clean = ctx.clean
    handles = ctx.handles
    legend = ctx.legend
    panel = ctx.panel
    save = ctx.save
    topology_legend = ctx.topology_legend
    os.environ['MPLCONFIGDIR'] = str(R / 'tmp/mplconfig')
    D = R / 'data/single_introduction'
    Q = R / 'qa'
    Q.mkdir(exist_ok=True)
    pairs = pd.read_csv(D / 'P5_PAIRED_OUTCOMES_192.csv')
    summary = json.loads((D / 'P5_SUMMARY_BOOTSTRAP_5000.json').read_text())
    rrg = pd.read_csv(D / 'P8_HIERARCHICAL_BOOTSTRAP_SUMMARY.csv').set_index('estimand')
    graphs = pd.read_csv(D / 'P8_FOUR_GRAPH_EFFECTS.csv')
    assert len(pairs) == 192 and len(graphs) == 4
    assert pairs[['i0A', 'i0F']].eq(1).all().all()
    rows = []
    means = []
    for beta in [20, 50]:
        p = pairs[pairs.beta == beta]
        assert len(p) == 96 and p.pair.nunique() == 96
        a, f = (p.rA.mean(), p.rF.mean())
        old = summary['observed'][f'AminusF_beta{beta}']
        assert np.isclose(a - f, old, atol=1e-14, rtol=0)
        l, u = summary['ci95'][f'beta{beta}']
        means.extend([dict(beta=beta, strategy='Adaptive', mean=float(a)), dict(beta=beta, strategy='Frozen', mean=float(f))])
        rows.append(dict(topology='4.8.8', beta=beta, estimate=float(f - a), lo=-u, hi=-l, old_estimate=old, old_lo=l, old_hi=u, graphs=1, pairs_per_graph=96, bootstrap=5000, interval_scope='fixed graph', source='P5_SUMMARY_BOOTSTRAP_5000.json'))
        v = rrg.loc[f'D{beta}']
        assert bool(v.CI_descriptive_exploratory)
        assert np.isclose(graphs[f'D{beta}'].mean(), v.point, atol=1e-14, rtol=0)
        rows.append(dict(topology='RRG3', beta=beta, estimate=-float(v.point), lo=-float(v.ci95_upper_graph_boot), hi=-float(v.ci95_lower_graph_boot), old_estimate=float(v.point), old_lo=float(v.ci95_lower_graph_boot), old_hi=float(v.ci95_upper_graph_boot), graphs=4, pairs_per_graph=16, bootstrap=10000, interval_scope='exploratory graph-first bootstrap', source='P8_HIERARCHICAL_BOOTSTRAP_SUMMARY.csv'))
    table = pd.DataFrame(rows)
    table.to_csv(D / 'FROZEN_MINUS_ADAPTIVE_ESTIMATES.csv', index=False, float_format='%.17g')
    pd.DataFrame(means).to_csv(D / 'PAIRED_488_MEANS.csv', index=False, float_format='%.17g')
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.5))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.21, top=0.91, wspace=0.32)
    purple = '#6F3FA0'
    grey = '#9E9E9E'
    xs = np.arange(2)
    width = 0.29
    means = pd.DataFrame(means)
    for strategy, off, face, hatch in [('Adaptive', -width / 2, purple, None), ('Frozen', width / 2, '#FFFFFF', '////')]:
        values = means[means.strategy == strategy].sort_values('beta')['mean'].to_numpy()
        axes[0].bar(xs + off, values, width=width, facecolor=face, edgecolor=purple, hatch=hatch, linewidth=0.85, label=strategy)
    axes[0].set(xticks=xs, xticklabels=['20', '50'], ylim=(0, 1), yticks=[0, 0.25, 0.5, 0.75, 1], xlabel='$\\beta_N$', ylabel='$R_\\infty$')
    axes[0].legend(loc='upper left', frameon=True, fancybox=False, framealpha=1, edgecolor='.7', fontsize=8)
    segments = []
    for topo, off, color, marker in [('4.8.8', -0.065, purple, '8'), ('RRG3', 0.065, grey, 'o')]:
        t = table[table.topology == topo].sort_values('beta')
        est = t.estimate.to_numpy()
        lo = t.lo.to_numpy()
        hi = t.hi.to_numpy()
        art = axes[1].errorbar(xs + off, est, yerr=[est - lo, hi - est], fmt=marker, linestyle='none', color=color, mfc=color, mec=color, markersize=5, capsize=3, elinewidth=0.9, capthick=0.8, label=topo)
        for (_, row), seg in zip(t.iterrows(), art.lines[2][0].get_segments()):
            assert np.allclose(seg[:, 1], [row.lo, row.hi], atol=1e-14, rtol=0)
            segments.append(dict(topology=topo, beta=int(row.beta), estimate=float(row.estimate), lo=float(seg[0, 1]), hi=float(seg[1, 1])))
    axes[1].axhline(0, color='.3', lw=0.65, ls='--', zorder=0)
    axes[1].set(xticks=xs, xticklabels=['20', '50'], ylim=(-0.25, 0.95), yticks=[-0.2, 0, 0.2, 0.4, 0.6, 0.8], xlabel='$\\beta_N$', ylabel='$\\Delta R_\\infty$')
    axes[1].legend(loc='upper left', frameon=True, fancybox=False, framealpha=1, edgecolor='.7', fontsize=8)
    for ax, letter in zip(axes, ['A', 'B']):
        ax.text(0, 1.045, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')
        ax.grid(axis='y', alpha=0.16, lw=0.4)
        ax.set_axisbelow(True)
        ax.set_xlim(-0.45, 1.45)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for txt in fig.findobj(Text):
        if not txt.get_visible() or not txt.get_text():
            continue
        b = txt.get_window_extent(renderer)
        if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > fig.bbox.width + 0.5 or (b.y1 > fig.bbox.height + 0.5):
            outside.append(txt.get_text())
    assert not outside, outside
    out = R / 'supplement/figures/SINGLE_INTRO_UPDATING_TWO_PANELS.pdf'
    ctx.export(fig, out, metadata={'CreationDate': None, 'ModDate': None, 'Creator': 'FPM two-panel paired-estimate renderer'})
    ctx.export(fig, Q / 'SINGLE_INTRO_PREVIEW.png', dpi=180)
    plt.close(fig)
    qa = dict(status='PASS', source_data_rows=192, graph_level_rows=4, new_simulations=0, new_bootstrap_replicates=0, source_means_replayed=True, source_intervals_reused=True, convention='Frozen - Adaptive; [lo,hi] = [-old_hi,-old_lo]', segments=segments, text_outside_page=outside, font='STIXGeneral, matching CAS body family', topology_style={'4.8.8': {'color': purple, 'marker': '8'}, 'RRG3': {'color': grey, 'marker': 'o'}}, topology_legend_order=['4.8.8', 'RRG3'], panels=2, caption_in_latex=True, raster_backgrounds=0, pdf_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
    (Q / 'SINGLE_INTRO_RENDER.json').write_text(json.dumps(qa, indent=2) + '\n', encoding='utf-8')
