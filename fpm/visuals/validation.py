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

from .style import ORDER as NETWORK_ORDER, STYLE as NETWORK_STYLE
import pymupdf as fitz
from matplotlib.colors import ListedColormap
from statistics import NormalDist

def retained_source_regression(root):
    R = Path(root) / 'reproduction/si'
    ORDER = NETWORK_ORDER
    STYLE = NETWORK_STYLE

    def vertices(path):
        items = path['items']
        assert all((v[0] == 'l' for v in items))
        return np.array([tuple(items[0][1])] + [tuple(v[2]) for v in items])
    report = {}
    page = fitz.open(R / 'reference/figures_before/SUPP_FIG_S4_ENDPOINT_THRESHOLD_ROBUSTNESS.pdf')[0]
    paths = [p for p in page.get_drawings() if len(p['items']) == 40]
    assert len(paths) == 9
    d = pd.read_csv(R / 'data/s1/ENDPOINT_THRESHOLD_VALUES.csv')
    errors = []
    for j, obs in enumerate(['T_epi', 'reach', 'q']):
        vals = []
        coords = []
        for k, factor in enumerate([0.5, 1, 2]):
            t = d[d.factor == factor].sort_values('beta_n')
            xy = vertices(paths[3 * j + k])
            assert len(xy) == 41
            fx = np.polyfit(t.beta_n, xy[:, 0], 1)
            errors.append(float(np.max(np.abs(np.polyval(fx, t.beta_n) - xy[:, 0]))))
            vals.extend(t[obs])
            coords.extend(xy[:, 1])
        fy = np.polyfit(vals, coords, 1)
        errors.append(float(np.max(np.abs(np.polyval(fy, vals) - coords))))
    assert max(errors) < 0.0002, max(errors)
    report['S1'] = dict(status='PASS', source_rows=123, original_vector_curves=9, plotted_coordinates=369, max_PDF_coordinate_residual_pt=max(errors), exact_source_mapping_closed=True)
    page = fitz.open(R / 'reference/figures_before/F05F_HEATMAP_ROBUST_PEAK_COUNT.pdf')[0]
    source = pd.read_csv(R / 'data/s2/PAPER_WINDOW_5OBS.csv')
    d = source[source.observable == 'robust_peak_count']
    cmap = ListedColormap(matplotlib.colormaps['viridis'](np.linspace(0.02, 0.98, 256)))
    cells = []
    for info, topo in zip(page.get_image_info(xrefs=True)[:8], ORDER):
        pix = fitz.Pixmap(page.parent, info['xref'])
        rgb = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3]
        z = d[d.topology == topo].pivot(index='delta', columns='beta_n', values='mean').sort_index().sort_index(axis=1).to_numpy()[::-1]
        ys = ((np.arange(20) + 0.5) * pix.height / 20).astype(int)
        xs = ((np.arange(28) + 0.5) * pix.width / 28).astype(int)
        actual = rgb[np.ix_(ys, xs)].astype(int)
        expected = cmap(z / 4, bytes=True)[:, :, :3].astype(int)
        err = int(np.abs(actual - expected).max())
        assert err <= 3, (topo, err)
        cells.append(dict(topology=topo, cells=560, max_RGB_channel_difference=err))
    report['S2'] = dict(status='PASS', measured_cells=4480, legacy_color_comparison='compatible at 8-bit palette quantization; not an exact numerical inversion', compatible_legacy_palette='viridis truncated to 0.02..0.98; inferred color compatibility, not recovered renderer metadata', new_palette='standard viridis, same 0..4 observable limits', panels=cells, interpolated_cells=0)
    d = pd.read_csv(R / 'data/s3/ORIGINAL_PAIRED_CONTRASTS.csv')
    d = d[d.metric == 'reach']
    a = pd.read_csv(R / 'data/s3/AUDITED_128_REACH_CONTRASTS.csv')
    t = d.merge(a, on=['topology', 'anchor', 'contrast'], suffixes=('', '_audit'), validate='one_to_one')
    assert len(t) == 128
    err = max((float(np.abs(t[c] - t[c + '_audit']).max()) for c in ['estimate', 'ci_low', 'ci_high']))
    assert err < 1e-12
    report['S3'] = dict(status='PASS', source_contrasts=128, displayed_contrasts=96, source_and_audited_intervals_max_difference=err, bootstrap_replicates_reused=1000, new_bootstrap_replicates=0)
    page = fitz.open(R / 'reference/figures_before/SUPP_FIG_I1_SINGLE_SEED_PMACRO_AND_DETECTOR.pdf')[0]
    d = pd.read_csv(R / 'data/s4/PRIMARY_832_COUNTS.csv').sort_values('beta_n')
    prob = (d.n_ge5pct / d.n).to_numpy()
    blue = lambda p: p['color'] is not None and abs(p['color'][0] - 0.1215686) < 0.002
    curves = [p for p in page.get_drawings() if blue(p) and len(p['items']) == 7]
    curve = next((p for p in curves if p['rect'].x0 < 100))
    xy = vertices(curve)
    fx = np.polyfit(np.log2(d.beta_n), xy[:, 0], 1)
    xerr = float(np.max(np.abs(np.polyval(fx, np.log2(d.beta_n)) - xy[:, 0])))
    assert xerr < 0.0002
    fy = np.polyfit(prob, xy[:, 1], 1)
    yerr = float(np.max(np.abs(np.polyval(fy, prob) - xy[:, 1])))
    assert yerr < 0.0002
    bars = [p for p in page.get_drawings() if blue(p) and len(p['items']) == 1 and (abs(p['rect'].width) < 0.0001) and (p['rect'].x0 < 420) and (p['rect'].height > 1)]
    assert len(bars) == 8
    z = NormalDist().inv_cdf(0.975)
    n = d.n.to_numpy()
    den = 1 + z * z / n
    center = (prob + z * z / (2 * n)) / den
    half = z * np.sqrt(prob * (1 - prob) / n + z * z / (4 * n * n)) / den
    bounds = np.array([[min((p['rect'].y0 - fy[1]) / fy[0], (p['rect'].y1 - fy[1]) / fy[0]), max((p['rect'].y0 - fy[1]) / fy[0], (p['rect'].y1 - fy[1]) / fy[0])] for p in bars])
    cierr = float(np.max(np.abs(bounds - np.column_stack([center - half, center + half]))))
    assert cierr < 1e-06
    report['S4'] = dict(status='PASS', total_realizations=int(d.n.sum()), probability_points=8, cutoff_comparison_points=24, original_axis='log2', selected_axis='linear; explicit author choice', log2_axis_PDF_residual_pt=xerr, probability_PDF_residual_pt=yerr, max_original_Wilson_bound_difference=cierr, conditional_Wilson_columns_used=False)
    report['status'] = 'PASS'
    (R / 'qa/SOURCE_REGRESSION.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('Source regression PASS: S1 curves, S2 measured cells, S3 original intervals, S4 counts and original Wilson bars.')
