"""Reproduce current measured-data figures through shared modules."""
from pathlib import Path
import hashlib,json,shutil
import pymupdf as fitz
from . import response,topology,sensitivity,article
from .style import PlotContext
from .inference import conditional_establishment
from .validation import retained_source_regression

JOBS=[
    ('state_scheme','response','standard','main','reproduction/main',response.state_scheme),
    ('epidemic_histories','response','standard','main','reproduction/main',response.epidemic_histories),
    ('peak_duration','response','standard','main','article_update',article.peak_duration),
    ('measured_maps','response','standard','main','reproduction/main',response.measured_maps),
    ('recurrence_depth','response','standard','main','reproduction/main',response.recurrence_depth),
    ('establishment_response','topology','establishment','main','reproduction/main/preserved_renderers/fig3',topology.establishment_response),
    ('population_crossovers','topology','crossovers','main','reproduction/main/preserved_renderers/fig7',topology.population_crossovers),
    ('response_landmarks','topology','landmarks','main','article_update',article.selected_landmarks),
    ('paired_controls','topology','standard','main','reproduction/main',topology.paired_controls),
    ('degree_matched_reach','topology','standard','main','article_update',article.degree_window),
    ('endpoint_sensitivity','sensitivity','supplement','si','reproduction/si',sensitivity.endpoint_sensitivity),
    ('mean_peak_map','sensitivity','supplement','si','reproduction/si',sensitivity.mean_peak_map),
    ('control_sensitivity','sensitivity','supplement','si','reproduction/si',sensitivity.control_sensitivity),
    ('single_introduction','sensitivity','supplement','si','reproduction/si',sensitivity.single_introduction),
    ('single_introduction_controls','sensitivity','paired','si','reproduction/si',sensitivity.single_introduction_controls),
    ('relative_transmission_scales','sensitivity','standard','si','article_update',article.relative_landmark_bars),
]

def render(root,group='all'):
    root=Path(root);(root/'qa').mkdir(exist_ok=True)
    index=root/'article_update/ACTIVE_FIGURE_ASSETS.json'
    if not index.is_file():raise FileNotFoundError('Attach the current figure update with prepare-data --update-archive before rendering.')
    expected=json.loads(index.read_text('utf-8'))
    if group in {'all','topology'}:conditional_establishment(root)
    commands=[];generated={'main':set(),'si':set()}
    for name,category,profile,namespace,relative,handler in JOBS:
        if group!='all' and group!=category:continue
        print(name,flush=True)
        with PlotContext(root/relative,profile,namespace) as context:
            handler(context)
            destination=root/'reproduction'/namespace/('manuscript/figures' if namespace=='main' else 'supplement/figures')
            destination.mkdir(parents=True,exist_ok=True)
            for p in context.exports:
                if p.suffix!='.pdf':continue
                target=destination/p.name
                if p.resolve()!=target.resolve():shutil.copyfile(p,target)
                generated[namespace].add(p.name)
        commands.append(dict(recipe=name,module=handler.__module__,group=category))
    if group in {'all','sensitivity'}:retained_source_regression(root)
    comparisons=[]
    for row in expected:
        namespace=row['document'];name=row['name']
        if name not in generated[namespace]:continue
        produced=root/row['generated'];reference=root/row['reference']
        assert produced.is_file(),row['generated']
        actual=hashlib.sha256(produced.read_bytes()).hexdigest()
        assert hashlib.sha256(reference.read_bytes()).hexdigest()==row['sha256']
        with fitz.open(reference) as a,fitz.open(produced) as b:
            pixels=len(a)==len(b) and all(x.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).samples==y.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).samples for x,y in zip(a,b))
        comparisons.append(dict(figure_asset=name,document=namespace,reference_sha256=row['sha256'],
            regenerated_sha256=actual,byte_identical=row['sha256']==actual,pixels_identical=pixels))
        assert pixels and row['sha256']==actual,name
    assert len(comparisons)=={'all':18,'response':7,'topology':5,'sensitivity':6}[group]
    receipt=dict(status='PASS_UNIFIED_VISUAL_MODULES',group=group,commands=commands,figure_assets=len(comparisons),
        comparisons=comparisons,legacy_renderer_execution=False,legacy_renderer_import=False,
        new_epidemics=0,new_ODE_campaign=0,paper_overwritten=False)
    (root/'qa/ALL_FIGURES_REBUILD.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'],figure_assets=len(comparisons),
        all_byte_identical=all(x['byte_identical'] for x in comparisons))))
    return receipt
