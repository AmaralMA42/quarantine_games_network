"""Render the approved measured-data figures through ordinary module calls."""
from pathlib import Path
import hashlib,json,shutil
import pymupdf as fitz
from . import response,topology,sensitivity
from .style import PlotContext
from .inference import conditional_establishment
from .validation import retained_source_regression

# Recipe names describe scientific content. The same dynamics implementation
# serves every experiment; none of these rendering jobs simulates an epidemic.
JOBS=[
    ('state_scheme','response','standard','main','reproduction/main',response.state_scheme),
    ('epidemic_histories','response','standard','main','reproduction/main',response.epidemic_histories),
    ('parameter_response','response','standard','main','reproduction/main',response.parameter_response),
    ('measured_maps','response','standard','main','reproduction/main',response.measured_maps),
    ('recurrence_depth','response','standard','main','reproduction/main',response.recurrence_depth),
    ('establishment_response','topology','establishment','main','reproduction/main/preserved_renderers/fig3',topology.establishment_response),
    ('population_crossovers','topology','crossovers','main','reproduction/main/preserved_renderers/fig7',topology.population_crossovers),
    ('response_landmarks','topology','landmarks','main','reproduction/main/selected_fig8',topology.response_landmarks),
    ('paired_controls','topology','standard','main','reproduction/main',topology.paired_controls),
    ('degree_matched_reach','topology','standard','main','reproduction/main',topology.degree_matched_reach),
    ('endpoint_sensitivity','sensitivity','supplement','si','reproduction/si',sensitivity.endpoint_sensitivity),
    ('mean_peak_map','sensitivity','supplement','si','reproduction/si',sensitivity.mean_peak_map),
    ('control_sensitivity','sensitivity','supplement','si','reproduction/si',sensitivity.control_sensitivity),
    ('single_introduction','sensitivity','supplement','si','reproduction/si',sensitivity.single_introduction),
    ('single_introduction_controls','sensitivity','paired','si','reproduction/si',sensitivity.single_introduction_controls),
]
RELOCATE={
    'establishment_response':('figures/FIG3_CONDITIONAL_AB_MEANS_CD.pdf','FIG3_CONDITIONAL_AB_MEANS_CD.pdf'),
    'population_crossovers':('figures/FIG7_SOLID_GUIDES_RRG_LAST.pdf','FIG7_SOLID_GUIDES_RRG_LAST.pdf'),
    'response_landmarks':('figures/FIG8_BROKEN_LINEAR_POINTS_ONLY.pdf','FIG8_SELECTED_BROKEN_LINEAR.pdf'),
}

def render(root,group='all'):
    root=Path(root);(root/'qa').mkdir(exist_ok=True)
    if group in {'all','topology'}:conditional_establishment(root)
    commands=[];generated={'main':set(),'si':set()}
    for name,category,profile,namespace,relative,handler in JOBS:
        if group!='all' and group!=category:continue
        print(name,flush=True)
        with PlotContext(root/relative,profile,namespace) as context:
            handler(context)
            generated[namespace].update(p.name for p in context.exports if p.suffix=='.pdf')
            if name in RELOCATE:
                source,target=RELOCATE[name]
                destination=root/'reproduction/main/manuscript/figures'/target
                destination.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(context.root/source,destination)
                generated['main'].add(target)
        commands.append(dict(recipe=name,module=handler.__module__,group=category))
    if group in {'all','sensitivity'}:retained_source_regression(root)
    comparisons=[]
    for folder,namespace in [('manuscript','main'),('supplement','si')]:
        if group!='all' and ((namespace=='si')!=(group=='sensitivity')):continue
        for approved in sorted((root/'paper'/folder/'figures').glob('*.pdf')):
            if approved.name not in generated[namespace]:continue
            produced=root/'reproduction'/namespace/folder/'figures'/approved.name
            if not produced.exists():continue
            expected=hashlib.sha256(approved.read_bytes()).hexdigest();actual=hashlib.sha256(produced.read_bytes()).hexdigest()
            a,b=fitz.open(approved),fitz.open(produced)
            pixels=len(a)==len(b) and all(x.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).samples==y.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).samples for x,y in zip(a,b))
            comparisons.append(dict(figure_asset=approved.name,document=folder,approved_sha256=expected,regenerated_sha256=actual,byte_identical=expected==actual,pixels_identical=pixels))
            assert pixels and expected==actual,approved.name
    assert len(comparisons)=={'all':17,'response':7,'topology':5,'sensitivity':5}[group]
    receipt=dict(status='PASS_UNIFIED_VISUAL_MODULES',group=group,commands=commands,figure_assets=len(comparisons),comparisons=comparisons,legacy_renderer_execution=False,legacy_renderer_import=False,new_epidemics=0,new_ODE_campaign=0,paper_overwritten=False)
    (root/'qa/ALL_FIGURES_REBUILD.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=receipt['status'],figure_assets=len(comparisons),all_byte_identical=all(x['byte_identical'] for x in comparisons))))
    return receipt
