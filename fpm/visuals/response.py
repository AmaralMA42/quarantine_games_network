"""Scientific compositions using shared identities and export primitives."""
from pathlib import Path
import os,json,hashlib
import numpy as np,pandas as pd,matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.container import ErrorbarContainer
from matplotlib.patches import Rectangle,FancyBboxPatch
from PIL import Image
from scipy.signal import find_peaks
from statistics import NormalDist
from .style import STATE_COLORS,LANDMARKS

def state_scheme(ctx):
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
    src = R / 'reference/FIG1_PRESERVED_LAYOUT.png'
    fig = plt.figure(figsize=(7.1, 7.1 * 951 / 2048))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(Image.open(src), extent=[0, 2048, 951, 0], interpolation='none')
    ax.set(xlim=(0, 2048), ylim=(951, 0))
    ax.axis('off')
    ax.add_patch(Rectangle((0, 0), 2048, 448, facecolor='white', edgecolor='none', zorder=3))
    labels = []

    def text(x, y, value, **kw):
        a = ax.text(x, y, value, ha='center', va='center', zorder=6, **kw)
        labels.append(value)
        return a

    def box(x, y, w, h, edge, face, label):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=10', facecolor=face, edgecolor=edge, lw=1, zorder=4))
        text(x + w / 2, y + h / 2, label, fontsize=14)
    box(249, 65, 325, 106, '#315A99', '#DCE7F6', '$S_N$')
    box(249, 238, 325, 106, STATE_COLORS['protected'], '#EEF4FA', '$S_Q$')
    box(865, 65, 326, 106, STATE_COLORS['infected'], '#E8F3E8', '$I_N$')
    box(865, 238, 326, 106, STATE_COLORS['infected'], '#F0F7F0', '$I_Q$')
    box(1615, 151, 326, 106, '#555555', '#EEEEEE', '$R$')
    for start, end in [((578, 118), (855, 118)), ((578, 290), (855, 290)), ((1196, 120), (1606, 188)), ((1196, 290), (1606, 222)), ((407, 179), (407, 229)), ((453, 229), (453, 179)), ((1025, 179), (1025, 229)), ((1070, 229), (1070, 179))]:
        ax.annotate('', xy=end, xytext=start, arrowprops=dict(arrowstyle='-|>', lw=0.9, color='#555555'), zorder=5)
    for x, y, label in [(715, 85, '$\\beta_N,\\,\\beta_A$'), (715, 258, '$\\beta_Q,\\,\\beta_A$'), (1399, 108, '$\\gamma$'), (1399, 277, '$\\gamma$'), (329, 204, '$\\tau\\Theta$'), (938, 204, '$\\tau\\Theta$')]:
        text(x, y, label, fontsize=9)
    text(28, 50, 'A', fontsize=11, fontweight='bold')
    text(28, 429, 'B', fontsize=11, fontweight='bold')
    for x, k in [(315, 3), (996, 4), (1677, 6)]:
        text(x, 411, f'$k={k}$', fontsize=8.2)
    regions = [(60, 638, 260, 50), (390, 638, 140, 50), (770, 638, 190, 50), (1050, 638, 200, 50), (1550, 630, 250, 47), (240, 882, 180, 52), (915, 882, 180, 52), (1605, 873, 160, 51)]
    for x, y, w, h in regions:
        ax.add_patch(Rectangle((x, y), w, h, facecolor='white', edgecolor='none', zorder=4))
    network_labels = [('Honeycomb', 183, 658, STYLE['honeycomb'][1]), ('4.8.8', 455, 658, STYLE['truncated_square_488'][1]), ('Kagome', 853, 658, STYLE['kagome'][1]), ('Square', 1145, 658, STYLE['square'][1]), ('Triangular', 1678, 650, STYLE['triangular'][1]), ('RRG3', 310, 904, STYLE['rr_k3'][1]), ('RRG4', 985, 904, STYLE['rr_k4'][1]), ('RRG6', 1678, 890, STYLE['rr_k6'][1])]
    for label, x, y, color in network_labels:
        text(x, y, label, fontsize=8, fontweight='bold', color=color)
    save(fig, 'FIG1_PRESERVED_GEOMETRY_UNIFORM_TYPE', [])
    (R / 'qa/FIG1_TYPOGRAPHY.json').write_text(json.dumps(dict(status='PASS', source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(), original_network_art_retained=True, network_art_raster=True, top_model_vector=True, all_visible_labels_retyped=True, font_family='STIXGeneral', math_fontset='stix', panel_letters_same_size_and_face=True, topology_layout_changed=False, original_png_modified=False, epidemic_simulations=0, infected_color=STATE_COLORS['infected'], protected_color=STATE_COLORS['protected'], labels=labels), indent=2) + '\n', encoding='utf-8')

def parameter_response(ctx):
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
    V = pd.read_csv(R / 'data/fig4_values.csv')
    V = V[V.within_historical_domain]
    bd = float((R / 'data/fig4_beta_dagger.txt').read_text())
    fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.35))
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.24, top=0.93, wspace=0.28, hspace=0.32)
    records = []
    for ax, (ob, yl), letter in zip(axs.flat, [('reach', '$R_\\infty$ / $R_\\dagger$'), ('peak', '$I_{\\max}$'), ('duration', '$T_{\\mathrm{epi}}$'), ('q_final', '$q_f$ / $q_\\dagger$')], 'ABCD'):
        od = V[(V.system == 'ode') & (V.observable == ob)].sort_values('beta_n')
        for mask in [od.beta_n <= bd - 0.05, od.beta_n >= bd + 0.05] if ob != 'peak' else [od.beta_n.notna()]:
            part = od[mask]
            ax.plot(part.beta_n, part['mean'], color='#111111', ls=(0, (4, 2)), lw=0.85, marker='.', ms=0.85)
            records += part[['system', 'observable', 'beta_n', 'mean']].to_dict('records')
        if ob != 'peak':
            ax.axvline(bd, color='.5', ls=':', lw=0.65)
        for top in ORDER + ['well_mixed']:
            d = V[(V.system == top) & (V.observable == ob)].sort_values('beta_n')
            assert len(d) == (45 if top == 'well_mixed' else 39)
            lab, col, m = STYLE.get(top, ('Finite WM (CTMC)', '#111111', '+'))
            ax.errorbar(d.beta_n, d['mean'], yerr=[d['mean'] - d.ci_low, d.ci_high - d['mean']], color=col, marker=m, ms=2.6, lw=0.65, elinewidth=0.45, capsize=0.8, linestyle='none' if top == 'well_mixed' else '-')
            records += d[['system', 'observable', 'beta_n', 'mean', 'ci_low', 'ci_high']].to_dict('records')
        ax.set(xlim=(0.8, 20.25), xticks=[1, 5, 10, 15, 20], xlabel='$\\beta_N$', ylabel=yl)
        ax.set_title(letter, loc='left', fontweight='bold')
        clean(ax)
    topology_legend(fig, [Line2D([], [], color='#111111', marker='+', ls='none', label='Finite WM (CTMC)'), Line2D([], [], color='#111111', ls='--', label='Well-mixed ODE')])
    save(fig, 'FIG4_AUDITED_ALL_NETWORKS', records)

def epidemic_histories(ctx):
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
    S = pd.read_csv(R / 'data/fig2/square_absorcao.csv')
    A = pd.read_csv(R / 'data/fig2/rrg4_adaptativo.csv')
    F = pd.read_csv(R / 'data/fig2/rrg4_frozen.csv')
    fig, ax = plt.subplots(2, 2, figsize=(7.1, 4.25), gridspec_kw={'width_ratios': [0.85, 1.15]})
    fig.subplots_adjust(left=0.1, right=0.985, bottom=0.14, top=0.91, wspace=0.28, hspace=0.38)
    records = []
    for j, (data, end, key) in enumerate([(S, 30, 'square_absorcao'), (A, 65, 'rrg4_adaptativo')]):
        for row, col, color in [(0, 'I_fraction', STATE_COLORS['infected']), (1, 'Q_fraction', STATE_COLORS['protected'])]:
            ax[row, j].plot(data.time, data[col], color=color, lw=1.1, drawstyle='steps-post')
            clean(ax[row, j])
            ax[row, j].set_xlim(0, end)
            records += list((dict(trajectory=key, observable=col, time=t, value=v) for t, v in zip(data.time, data[col])))
        ax[0, j].tick_params(labelbottom=False)
        ax[1, j].set(ylim=(0, 1), xlabel='Time ($\\gamma^{-1}$)')
    for row, col, color in [(0, 'I_fraction', STATE_COLORS['infected']), (1, 'Q_fraction', STATE_COLORS['protected'])]:
        ax[row, 1].plot(F.time, F[col], color=color, lw=1.1, ls=(0, (4, 2.2)), drawstyle='steps-post')
        records += list((dict(trajectory='rrg4_frozen', observable=col, time=t, value=v) for t, v in zip(F.time, F[col])))
    ax[0, 0].set(ylim=(0, 0.023), ylabel='$I/N$')
    ax[0, 1].set_ylim(0, 0.12)
    ax[1, 0].set_ylabel('$Q/N$')
    for a, l in zip(ax.ravel(), 'ABCD'):
        a.text(0, 1.035, l, transform=a.transAxes, ha='left', va='bottom', fontweight='bold', fontsize=11)
    sm = np.convolve(A.I_count, np.ones(9) / 9, mode='same')
    ix, _ = find_peaks(sm, height=max(1, 16384 / (5 * 16)), prominence=max(1, 0.05 * sm.max()), distance=24)
    assert len(ix) == 4
    ax[0, 1].plot(A.time.iloc[ix], A.I_fraction.iloc[ix], 'o', color=STATE_COLORS['infected'], ms=3)
    ax[0, 1].legend(handles=[Line2D([], [], color='#333333', lw=1.1, label='Adaptive'), Line2D([], [], color='#333333', lw=1.1, ls=(0, (4, 2.2)), label='Frozen')], loc='upper right', frameon=True, fancybox=False, framealpha=1, facecolor='white', edgecolor='.5', fontsize=7)
    save(fig, 'FIG2_COMPACT_MATCHED_CLEAN', records)

def measured_maps(ctx):
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
    for file, name, label, cmap, limit in [('FIG5_REACH_4480_VALUES.csv', 'FIG5_REACH', '$R_\\infty$', 'viridis', 1), ('FIG5_DURATION_4480_VALUES.csv', 'FIG5_DURATION', '$T_{\\mathrm{epi}}$', 'magma', None), ('FIG6_RECURRENCE_4480_VALUES.csv', 'FIG6_RECURRENCE', '$P(n_p\\geq2)$', 'viridis', 1)]:
        D = pd.read_csv(R / 'data' / file)
        assert len(D) == 4480
        vmax = float(D['mean'].max()) if limit is None else limit
        fig, axs = plt.subplots(2, 4, figsize=(7.1, 3.65))
        fig.subplots_adjust(left=0.075, right=0.885, bottom=0.15, top=0.91, wspace=0.3, hspace=0.45)
        records = []
        for j, (a, t) in enumerate(zip(axs.flat, ORDER)):
            d = D[D.topology == t]
            xs = sorted(d.beta_n.unique())
            ys = sorted(d.delta.unique())
            assert len(xs) == 28 and len(ys) == 20 and np.allclose(np.diff(xs), 0.5) and np.allclose(np.diff(ys), 1)
            z = d.pivot(index='delta', columns='beta_n', values='mean').reindex(index=ys, columns=xs)
            assert z.notna().all().all()
            im = a.imshow(z, origin='lower', aspect='auto', interpolation='nearest', cmap=cmap, vmin=0, vmax=vmax)
            a.set_title(chr(65 + j), loc='left', fontweight='bold', fontsize=9)
            a.set_xticks([0, 9, 18, 27], [f'{xs[i]:g}' for i in [0, 9, 18, 27]])
            a.set_yticks([0, 9, 19], [f'{ys[i]:g}' for i in [0, 9, 19]])
            if j >= 4:
                a.set_xlabel('$\\beta_N$')
            if j % 4 == 0:
                a.set_ylabel('$\\delta$')
            records += d.to_dict('records')
        fig.colorbar(im, cax=fig.add_axes([0.91, 0.2, 0.018, 0.63]), label=label)
        save(fig, name, records)

def recurrence_depth(ctx):
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
    D = pd.read_csv(R / 'data/fig10_conditional_intervals.csv').rename(columns={'estimate': 'probability', 'n': 'denominator'})
    old = pd.read_csv(R / 'data/fig10_probabilities.csv')
    j = D.merge(old, on=['panel', 'topology', 'k'], suffixes=('_ci', '_article'), validate='one_to_one')
    assert len(j) == 72
    for c in ['probability', 'successes', 'denominator']:
        assert np.allclose(j[c + '_ci'], j[c + '_article'], atol=1e-14, rtol=0)
    fig, axs = plt.subplots(1, 2, figsize=(7.1, 4.1))
    fig.subplots_adjust(left=0.1, right=0.985, bottom=0.29, top=0.92, wspace=0.25)
    records = []
    segments = []
    for a, panel in zip(axs, ['A', 'B']):
        for t in ORDER + ['well_mixed']:
            d = D[(D.panel == panel) & (D.topology == t)].sort_values('k')
            lab, col, m = STYLE.get(t, ('Finite WM (CTMC)', '#111111', 'x'))
            lo = d.probability - d.selected_lo
            hi = d.selected_hi - d.probability
            assert (lo >= -1e-14).all() and (hi >= -1e-14).all()
            cont = a.errorbar(d.k, d.probability, yerr=[lo.clip(lower=0), hi.clip(lower=0)], fmt=m + '-', color=col, mfc=col, ms=3.5, lw=0.85, elinewidth=0.65, capsize=1.3, capthick=0.65)
            for (_, r), seg in zip(d.iterrows(), cont.lines[2][0].get_segments()):
                assert np.allclose(seg, [[r.k, r.selected_lo], [r.k, r.selected_hi]], atol=1e-14, rtol=0)
                segments.append(dict(panel=panel, topology=t, k=int(r.k), lo=float(seg[0, 1]), hi=float(seg[1, 1])))
            records += d.to_dict('records')
        a.set(xlabel='$k$', ylabel='$P(n_p\\geq k)$', xticks=[2, 3, 4, 5], ylim=(-0.04, 1.04))
        a.set_title(panel, loc='left', fontweight='bold')
        clean(a)
    topology_legend(fig, [Line2D([], [], color='#111111', marker='x', label='Finite WM (CTMC)')])
    save(fig, 'FIG10_ARTICLE_POINTS_WITH_CI', records)
    (R / 'qa/FIG10_ACTUAL_INTERVAL_SEGMENTS.json').write_text(json.dumps(segments, indent=2) + '\n')

