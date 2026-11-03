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

def establishment_response(ctx):
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
    O.mkdir(parents=True, exist_ok=True)
    for row in json.loads((R / 'provenance/INPUT_BINDINGS.json').read_text()):
        assert hashlib.sha256((R / row['stored_path']).read_bytes()).hexdigest() == row['sha256']
    SP = ORDER[:5]
    RR = ORDER[5:]
    all_cells = pd.read_csv(R / 'data/FIG3_ALL_392_MEASURED_CELLS.csv')
    est = pd.read_csv(R / 'data/FIG3_ESTABLISHMENT_CONDITIONAL_392.csv')
    reach = pd.read_csv(R / 'data/FIG3_REACH_392_MEASURED_POINTS.csv')
    qfinal = pd.read_csv(R / 'data/FIG3_Q_FINAL_392_MEASURED_POINTS.csv')
    assert len(all_cells) == len(est) == len(reach) == len(qfinal) == 392
    assert not est.duplicated(['topology', 'beta_n']).any()
    assert set(est.campaign) == {'M27', 'AB'}
    assert len(est[est.campaign == 'AB']) == 40
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.0))
    fig.subplots_adjust(left=0.087, right=0.985, bottom=0.105, top=0.925, wspace=0.28, hspace=0.36)
    ledger = []
    guides = []

    def common(ax, letter, ylabel, xlim, xticks):
        ax.set(xlabel='$\\beta_N$', ylabel=ylabel, xlim=xlim, ylim=(-0.04, 1.045))
        ax.set_xticks(xticks)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        ax.grid(axis='y', alpha=0.13, lw=0.45)
        ax.set_axisbelow(True)
        ax.text(0, 1.035, letter, transform=ax.transAxes, ha='left', va='bottom', fontsize=12, fontweight='bold')
        ax.xaxis.labelpad = 3
        ax.yaxis.labelpad = 3

    def record(frame, panel, observable):
        show = panel in ['A', 'B']
        for row in frame.to_dict('records'):
            ledger.append(dict(panel=panel, observable=observable, cell_id=row['cell_id'], campaign=row['campaign'], wave=row['wave'], topology=row['topology'], beta_n=row['beta_n'], estimate=row['estimate'], n=row['n'], n_graphs=row['n_graphs'], dynamics_per_graph=row['dynamics_per_graph'], marker_fill='filled', source_ci_lo=row['source_ci_lo'], source_ci_hi=row['source_ci_hi'], source_ci_method=row['source_ci_method'], ci_drawn=show, ci_drawn_lo=row['conditional_lo'] if show else None, ci_drawn_hi=row['conditional_hi'] if show else None, selected_ci_method=row['conditional_method'] if show else 'NOT_DISPLAYED_BY_HUMAN_REQUEST', ci_zero_width=False if show else None))
    artist_intervals = []

    def draw_intervals(ax, frame, color, panel):
        lo = frame.estimate - frame.conditional_lo
        hi = frame.conditional_hi - frame.estimate
        assert np.isfinite(np.vstack([lo, hi])).all()
        assert (lo >= -1e-14).all() and (hi >= -1e-14).all()
        c = ax.errorbar(frame.beta_n, frame.estimate, yerr=np.vstack([lo.clip(lower=0), hi.clip(lower=0)]), fmt='none', ecolor=color, elinewidth=0.65, capsize=1.6, capthick=0.65, zorder=6)
        segments = c.lines[2][0].get_segments()
        assert len(segments) == len(frame)
        for (_, row), seg in zip(frame.iterrows(), segments):
            assert np.allclose(seg, [[row.beta_n, row.conditional_lo], [row.beta_n, row.conditional_hi]], atol=1e-14, rtol=0)
            artist_intervals.append(dict(panel=panel, cell_id=row.cell_id, x_lo=float(seg[0, 0]), x_hi=float(seg[1, 0]), ci_lo=float(seg[0, 1]), ci_hi=float(seg[1, 1])))
    for ax, keys, panel in [(axes[0, 0], SP, 'A'), (axes[0, 1], RR, 'B')]:
        for topology in keys:
            label, color, marker = STYLE[topology]
            frame = est[(est.topology == topology) & (est.beta_n <= 2.6)].sort_values('beta_n')
            assert len(frame) == 24 and sum(frame.campaign == 'M27') == 19 and (sum(frame.campaign == 'AB') == 5)
            ax.plot(frame.beta_n, frame.estimate, color=color, lw=1.05, zorder=2)
            for campaign in ['M27', 'AB']:
                sub = frame[frame.campaign == campaign]
                ax.plot(sub.beta_n, sub.estimate, ls='none', marker=marker, ms=3.0 if campaign == 'M27' else 4.0, mfc=color, mec=color, mew=0.4, zorder=4)
            draw_intervals(ax, frame, color, panel)
            record(frame, panel, 'p_macro')
            guides.append(dict(panel=panel, topology=topology, coordinates=frame.cell_id.tolist(), beta=frame.beta_n.tolist(), y=frame.estimate.tolist()))
        common(ax, panel, '$P_{\\mathrm{macro}}$', (1, 2.6), [1, 1.4, 1.8, 2.2, 2.6])
        handles = [Line2D([], [], color=STYLE[t][1], marker=STYLE[t][2], mfc=STYLE[t][1], mec=STYLE[t][1], ms=3.7, lw=1, label=STYLE[t][0]) for t in (['truncated_square_488', 'honeycomb', 'kagome', 'square', 'triangular'] if panel == 'A' else keys)]
        ax.legend(handles=handles, loc='upper left' if panel == 'A' else 'lower right', frameon=False, fontsize=8.2, handlelength=1.35, handletextpad=0.45, labelspacing=0.23, borderpad=0.2)
    for ax, data, panel, obs, ylabel in [(axes[1, 0], reach, 'C', 'reach', '$R_\\infty$'), (axes[1, 1], qfinal, 'D', 'q_final', '$q_f$')]:
        for topology in ORDER:
            label, color, marker = STYLE[topology]
            frame = data[(data.topology == topology) & (data.campaign == 'M27') & (data.beta_n <= 30)].sort_values('beta_n')
            assert len(frame) == 42 and frame.beta_n.max() == 28
            ax.plot(frame.beta_n, frame.estimate, color=color, lw=0.95, marker=marker, ms=2.7, mfc=color, mec=color, mew=0.35, zorder=2)
            record(frame, panel, obs)
            guides.append(dict(panel=panel, topology=topology, coordinates=frame.cell_id.tolist(), beta=frame.beta_n.tolist(), y=frame.estimate.tolist()))
        common(ax, panel, ylabel, (1, 30), [1, 5, 10, 15, 20, 25, 30])
    late_errorbar_containers = {panel: len(ax.containers) for panel, ax in [('C', axes[1, 0]), ('D', axes[1, 1])]}
    assert late_errorbar_containers == {'C': 0, 'D': 0}
    name = 'FIG3_CONDITIONAL_AB_MEANS_CD'
    ctx.export(fig, O / (name + '.pdf'), metadata={'Title': 'Fig3 - measured establishment and epidemic response', 'Creator': 'FPM audited CSV renderer; candidate', 'CreationDate': None, 'ModDate': None})
    ctx.export(fig, O / (name + '.svg'), metadata={'Title': 'Fig3 - measured establishment and epidemic response', 'Creator': 'FPM audited CSV renderer; candidate', 'Date': None})
    ctx.export(fig, O / (name + '.png'), dpi=600)
    (O / 'ACTUAL_ERRORBAR_SEGMENTS.json').write_text(json.dumps(artist_intervals, indent=2) + '\n', encoding='utf-8')
    plt.close(fig)
    pd.DataFrame(ledger).to_csv(O / 'PLOTTED_COORDINATES.csv', index=False, float_format='%.17g')
    (O / 'GUIDE_VERTEX_LEDGER.json').write_text(json.dumps(guides, indent=2) + '\n', encoding='utf-8')
    receipt = dict(status='CANDIDATE_UNPROMOTED', panel_count=4, panel_map={'A': 'old A+C: spatial establishment', 'B': 'old B+D: RRG establishment', 'C': 'old E: cumulative reach', 'D': 'old F: final quarantine'}, physical_inches=[7.1, 5.0], axis_limits={'A': [1, 2.6], 'B': [1, 2.6], 'C': [1, 30], 'D': [1, 30]}, plotted_markers=len(ledger), establishment_M27_cells=152, establishment_AB_cells=40, reach_M27_cells=336, q_final_M27_cells=336, late_last_measured_beta_displayed=28, measured_beta_30_created=False, late_omitted_coordinates=[32, 40], all_source_cells_retained=392, establishment_one_sorted_guide_per_topology=True, samples_pooled=False, smoothing=False, fitting=False, new_simulations=0, PNG_inputs=0, mission_subtitles_inside_figure=False, figure_footer=False, global_topology_legend=False, topology_legends_inside_upper_panels=True, confidence_intervals='A/B only: 120 exact CP spatial and 72 conditional heterogeneous Bernoulli Chernoff-KL observed RRG bank; C/D point estimates only; source intervals retained.', late_errorbar_containers=late_errorbar_containers, style=STYLE)
    (O / 'RENDER_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')

def population_crossovers(ctx):
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
    O = R / 'figures'
    P = R / 'qa/figure_previews'
    P.mkdir(parents=True, exist_ok=True)
    A = pd.read_csv(R / 'data/fig7_windows.csv')
    for kind, file in [('SOLID_RRG_LAST', 'fig7_archived_widths.csv')]:
        B = pd.read_csv(R / 'data' / file)
        fig, axs = plt.subplots(1, 2, figsize=(7.1, 4.0), gridspec_kw={'width_ratios': [0.95, 1.15]})
        fig.subplots_adjust(left=0.135, right=0.985, bottom=0.23, top=0.94, wspace=0.32)
        records = []
        for j, top in enumerate(ORDER):
            lab, col, m = STYLE[top]
            d = A[A.topology == top].set_index('landmark')
            axs[0].plot([d.loc['beta05', 'estimate'], d.loc['beta95', 'estimate']], [j, j], color=col, lw=2.5, alpha=0.7)
            axs[0].plot(d.loc['beta50', 'estimate'], j, marker=m, color=col, ms=4, ls='none')
            records += A[(A.topology == top) & A.landmark.isin(['beta05', 'beta50', 'beta95'])].to_dict('records')
            if kind == 'SOLID_RRG_LAST':
                last = B[(B.topology == top) & (B.campaign == 'M25_P2')].sort_values('population').iloc[-1]
                newest = B[(B.topology == top) & (B.campaign == 'M27_G3')].iloc[0]
                axs[1].plot(np.log2([last.population, newest.population]), [last.estimate, newest.estimate], color=col, lw=0.7, ls='-', zorder=1)
            for campaign, g in B[B.topology == top].groupby('campaign'):
                g = g.sort_values('population')
                xx = np.log2(g.population)
                if campaign == 'M25_P2':
                    axs[1].plot(xx, g.estimate, color=col, lw=0.7)
                for x, (_, r) in zip(xx, g.iterrows()):
                    conditional = r.valid_fraction < 0.95
                    axs[1].plot([x, x], [r.ci_low, r.ci_high], color=col, lw=0.8, ls=':' if conditional else '-')
                    axs[1].plot([x, x], [r.ci_low, r.ci_high], ls='none', marker='_', markersize=4, color=col, markeredgewidth=0.75)
                    axs[1].plot(x, r.estimate, marker=m, color=col, ms=4, ls='none', mfc=col)
                    if conditional:
                        axs[1].annotate('$\\dagger$', (x, r.estimate), xytext=(4, 3), textcoords='offset points', fontsize=7, color=col)
                records += g.to_dict('records')
        axs[0].set_yticks(range(8), [STYLE[t][0] for t in ORDER])
        axs[0].invert_yaxis()
        axs[0].set(xlim=(1.03, 2.4), xlabel='$\\beta_N$')
        axs[0].set_title('A', loc='left', fontweight='bold')
        axs[1].set(xlim=(13.8, 18.25), ylim=(0, 0.71), xlabel='Population $N$ (log$_2$ axis)', ylabel='$W_\\beta$')
        axs[1].set_xticks(range(14, 19), ['$2^{%d}$' % i for i in range(14, 19)])
        axs[1].set_title('B', loc='left', fontweight='bold')
        for ax in axs:
            clean(ax)
        legend(fig, handles(ORDER[:5]), 5, loc='lower center', bbox_to_anchor=(0.52, 0.075))
        legend(fig, handles(ORDER[5:]), 3, loc='lower center', bbox_to_anchor=(0.52, 0.015))
        save(fig, 'FIG7_SOLID_GUIDES_RRG_LAST', records)
    (R / 'qa/PLOTTED_COORDINATES.json').write_text(json.dumps(coordinates, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (R / 'qa/RENDER_QA.json').write_text(json.dumps(qa, indent=2) + '\n', encoding='utf-8')

def response_landmarks(ctx):
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
    data_path = ROOT / 'data/fig8_landmarks.csv'
    df = pd.read_csv(data_path)
    assert len(df) == 48 and (not df.duplicated(['topology', 'landmark']).any())
    assert (df.beta_n.between(1, 2.5) | df.beta_n.between(8, 35)).all()
    fig, (left, right) = plt.subplots(1, 2, sharey=True, figsize=(7.1, 4.0), gridspec_kw={'width_ratios': [1, 4.25]})
    fig.subplots_adjust(left=0.155, right=0.985, bottom=0.14, top=0.975, wspace=0.07)
    left.set_xlim(1, 2.5)
    right.set_xlim(8, 35)
    left.set_xticks([1, 1.5, 2, 2.5], ['1', '1.5', '2', '2.5'])
    right.set_xticks([8, 10, 15, 20, 25, 30, 35])
    left.set_yticks(range(8), NAMES)
    right.tick_params(axis='y', which='both', left=False, labelleft=False)
    right.spines['left'].set_visible(False)
    for ax in (left, right):
        ax.set_ylim(7.55, -0.55)
        ax.grid(axis='x', alpha=0.13, lw=0.4)
        ax.set_axisbelow(True)
    ledger = []
    for k, (key, (label, color, marker)) in enumerate(LANDMARKS.items()):
        for j, topology in enumerate(ORDER):
            row = df[(df.topology == topology) & (df.landmark == key)].iloc[0]
            ax = left if row.beta_n <= 2.5 else right
            y = j + (k - 2.5) * 0.1
            size = 4.4 if marker == '*' else 3.0
            artist = ax.plot([row.beta_n], [y], marker=marker, color=color, linestyle='none', markersize=size, markeredgewidth=0.65, markerfacecolor='white' if key == 'rebound_ci_separated' else color)[0]
            assert float(artist.get_xdata()[0]) == float(row.beta_n)
            ledger.append({'topology': topology, 'landmark': key, 'beta_n': float(row.beta_n), 'y': y, 'axis_segment': '1_to_2.5' if ax is left else '8_to_35', 'marker': marker, 'color': color, 'errorbar_drawn': False})
    break_lines = [left.axvline(2.5, color='#666666', linestyle=':', linewidth=0.8, zorder=2), right.axvline(8, color='#666666', linestyle=':', linewidth=0.8, zorder=2)]
    group_lines = [ax.axhline(4.5, color='#666666', linestyle=':', linewidth=0.7, zorder=1.5) for ax in (left, right)]
    for ax, x in ((left, 1), (right, 0)):
        for y in (0, 1):
            ax.plot([x, x], [y, y], transform=ax.transAxes, marker=[(-1, -0.65), (1, 0.65)], markersize=7, markeredgewidth=0.8, color='#222222', linestyle='none', clip_on=False)
    items = [Line2D([], [], color=c, marker=m, linestyle='none', markersize=6.1 if m == '*' else 4.4, label=label, markerfacecolor='white' if key == 'rebound_ci_separated' else c) for key, (label, c, m) in LANDMARKS.items()]
    legend = right.legend(handles=items, ncol=1, loc='lower right', frameon=True, fancybox=False, framealpha=1, facecolor='white', edgecolor='#777777', fontsize=7, handletextpad=0.5, labelspacing=0.4, borderpad=0.65, borderaxespad=0.7)
    legend.get_frame().set_linewidth(0.6)
    fig.supxlabel('$\\beta_N$', x=0.57, y=0.025, fontsize=9)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    box = legend.get_window_extent(renderer)
    axes_box = right.get_window_extent(renderer)
    assert axes_box.contains(box.x0, box.y0) and axes_box.contains(box.x1, box.y1)
    for row in ledger:
        if row['axis_segment'] == '8_to_35':
            x, y = right.transData.transform((row['beta_n'], row['y']))
            assert not box.expanded(1.02, 1.15).contains(x, y)
    outside = []
    for text in fig.findobj(matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            b = text.get_window_extent(renderer)
            if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > fig.bbox.width + 0.5 or (b.y1 > fig.bbox.height + 0.5):
                outside.append(text.get_text())
    assert not outside, outside
    assert not any((isinstance(c, ErrorbarContainer) for ax in (left, right) for c in ax.containers))
    assert all((line.get_linestyle() == 'None' for ax in (left, right) for line in ax.lines if line not in break_lines + group_lines))
    assert all((line.get_linestyle() == ':' for line in break_lines))
    ROOT.joinpath('figures').mkdir(exist_ok=True)
    ctx.export(fig, ROOT / 'figures/FIG8_BROKEN_LINEAR_POINTS_ONLY.pdf', metadata={'CreationDate': None, 'ModDate': None, 'Creator': 'FPM numeric candidate renderer'})
    ctx.export(fig, ROOT / 'figures/FIG8_BROKEN_LINEAR_POINTS_ONLY.png', dpi=200)
    plt.close(fig)
    pd.DataFrame(ledger).to_csv(ROOT / 'PLOTTED_COORDINATES.csv', index=False)
    qa = {'status': 'PASS', 'source_sha256': hashlib.sha256(data_path.read_bytes()).hexdigest(), 'source_records': 48, 'displayed_records': len(ledger), 'hidden_records': 0, 'source_intervals_retained': int(df.ci_low.notna().sum()), 'displayed_errorbars': 0, 'x_segments': [[1, 2.5], [8, 35]], 'scales': ['linear', 'linear'], 'unequal_segment_scale_explicit': True, 'reach_minimum_marker': '*', 'legend_inside_right_axes': True, 'legend_frame': True, 'legend_position': 'lower right within existing data space', 'x_label': 'beta_N only', 'dotted_vertical_break_boundaries': [2.5, 8], 'group_separator_y': 4.5, 'landmark_y_offset_step': 0.1, 'y_limits': [7.55, -0.55], 'figure_size_inches': [7.1, 4.0], 'marker_legend_overlap': False, 'text_outside_canvas': outside, 'epidemic_simulations': 0, 'main_manuscript_replaced': False}
    ROOT.joinpath('QA.json').write_text(json.dumps(qa, indent=2) + '\n', encoding='utf-8')

def paired_controls(ctx):
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
    E = pd.read_csv(R / 'data/fig9_contrasts.csv')
    fig = plt.figure(figsize=(7.1, 5.15))
    gs = fig.add_gridspec(2, 3, left=0.13, right=0.985, bottom=0.16, top=0.95, hspace=0.42, wspace=0.18)
    axs = [fig.add_subplot(gs[0, i]) for i in range(3)]
    records = []
    for ax, panel, title in zip(axs, ['A1', 'A2', 'A3'], ['A1  $\\beta_N=8$', 'A2  Reference minimum', 'A3  $\\beta_N=40$']):
        for j, top in enumerate(ORDER):
            r = E[(E.topology == top) & (E.main_panel == panel)].iloc[0]
            lab, col, m = STYLE[top]
            ax.errorbar(r.estimate, j, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]], fmt=m, color=col, ms=4, capsize=1)
            records.append(r.to_dict())
        ax.axvline(0, color='.3', ls='--', lw=0.7)
        ax.set_title(panel, loc='left', fontsize=9, fontweight='bold')
        ax.set_xlabel('$\\Delta R_\\infty$')
        ax.set_yticks(range(8), [STYLE[t][0] for t in ORDER] if ax == axs[0] else [])
        ax.invert_yaxis()
        clean(ax)
    ax = fig.add_subplot(gs[1, :])
    for j, top in enumerate(ORDER):
        for off, panel in [(-0.13, 'B_RANGE'), (0.13, 'B_ALIGNMENT')]:
            r = E[(E.topology == top) & (E.main_panel == panel)].iloc[0]
            lab, col, m = STYLE[top]
            ax.errorbar(r.estimate, j + off, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]], fmt=m, color=col, mfc=col if off < 0 else 'white', ms=4, capsize=1)
            records.append(r.to_dict())
    ax.axvline(0, color='.3', ls='--', lw=0.7)
    ax.set_yticks(range(8), [STYLE[t][0] for t in ORDER])
    ax.invert_yaxis()
    ax.set_title('B', loc='left', fontsize=9, fontweight='bold')
    ax.set_xlabel('$\\Delta R_\\infty$')
    clean(ax)
    legend(fig, [Line2D([], [], marker='o', color='.2', ls='none', label='Local - Global'), Line2D([], [], marker='o', color='.2', mfc='white', ls='none', label='Local - Shuffled')], 2, loc='lower center', bbox_to_anchor=(0.55, 0.005))
    save(fig, 'FIG9_AUDITED_CONTRASTS', records)

def degree_matched_reach(ctx):
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
    C = pd.read_csv(R / 'data/fig11_contrasts.csv')
    C = C[C.beta_n <= 20].copy()
    fig, ax = plt.subplots(figsize=(7.1, 3.8))
    fig.subplots_adjust(top=0.8, left=0.1, right=0.985, bottom=0.15)
    items = []
    records = []
    for top in ORDER[:5]:
        g = C[C.spatial == top].sort_values('beta_n')
        lab, col, m = STYLE[top]
        rr = STYLE[g.rrg.iloc[0]][0]
        ax.plot(g.beta_n, g.difference, color=col, marker=m, mfc=col, ms=3, lw=0.9, zorder=3)
        items.append(Line2D([], [], color=col, marker=m, label=lab + ' - ' + rr, ms=4))
        records += g.to_dict('records')
    ax.axhline(0, color='.2', ls='--', lw=0.7)
    ax.set(xlabel='$\\beta_N$', ylabel='$\\Delta^{\\mathrm{geom}}R_\\infty$', xlim=(0.5, 20.4), xticks=[1, 5, 10, 15, 20], ylim=(-0.55, 0.3))
    clean(ax)
    legend(fig, items, 3, loc='upper center')
    save(fig, 'FIG11_AUDITED_FIVE_PAIRS', records)
    assert len(records) == 195 and len(ax.containers) == 0 and (len(ax.collections) == 0)
    (R / 'qa/FIG11_MEANS_ONLY_RENDER.json').write_text(json.dumps(dict(status='PASS', plotted_estimates=195, plotted_curves=5, errorbar_containers=0, errorbar_collections=0, source_interval_columns_preserved=True, point_estimates_only=True, x_measured_max=20, y_limits=[-0.55, 0.3]), indent=2) + '\n')
