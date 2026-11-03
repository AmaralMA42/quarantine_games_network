"""Parameter sweeps are data consumed by one dynamics implementation."""
from pathlib import Path
import csv,itertools,json,os
from . import model

def make_sweep(template,configuration,output):
    template=Path(template).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    with template.open(newline='',encoding='utf-8') as f:bases=list(csv.DictReader(f,delimiter='\t'))
    settings=json.loads(Path(configuration).read_text(encoding='utf-8'))
    grid=settings['parameters'];repetitions=int(settings['repetitions'])
    if repetitions<1:raise ValueError('Positive repetitions required')
    allowed={'beta_n','beta_q','beta_mix_a','delta','omega','tau','fermi_k','q0','i0','behavior_mode','information_mode','information_radius','information_alpha'}
    if not set(grid)<=allowed:raise ValueError('Derived identity/graph fields cannot be sweep axes')
    if any(not isinstance(v,list) or not v for v in grid.values()):raise ValueError('Every parameter axis must be a nonempty list')
    paired=set(settings.get('paired_parameters',[]))
    if not paired<=set(grid):raise ValueError('Paired parameter must be a varied axis')
    rows=[];output.parent.mkdir(parents=True,exist_ok=True)
    for base_index,base in enumerate(bases):
        for values in itertools.product(*grid.values()):
            parameters=dict(zip(grid,values))
            for repeat in range(repetitions):
                row=dict(base);row.update({k:str(v) for k,v in parameters.items()})
                row['beta_a']=str(float(row['beta_q'])+float(row['beta_mix_a'])*(float(row['beta_n'])-float(row['beta_q'])))
                row['i_c']=str(float(row['omega'])/(float(row['delta'])*float(row['beta_n'])))
                identity=dict(master=settings['master_seed'],base=base_index,replicate=repeat,parameters={k:v for k,v in parameters.items() if k not in paired})
                master=model.stable_seed(json.dumps(identity,sort_keys=True,separators=(',',':')))
                for key,domain in [('q_init_seed','initial_Q'),('i_init_seed','initial_I'),('dynamic_seed','dynamics'),('signal_seed','signal')]:row[key]=str(model.derive_seed(master,domain))
                if row['topology']!='well_mixed':
                    graph=Path(base['graph_path'])
                    if not graph.is_absolute():graph=template.parent/graph
                    row['graph_path']=os.path.relpath(graph.resolve(),output.parent).replace(os.sep,'/')
                    row['graph_certificate_path']=row['graph_path']+'.periodic.json' if base['graph_certificate_path']!='-' else '-'
                row.update(run_index=str(len(rows)),run_id=f'parameter_sweep_{len(rows):06d}',schema_version=model.SCHEMA_VERSION,fingerprint_version=model.FINGERPRINT_VERSION)
                row['run_fingerprint']=model.semantic_run_fingerprint(row)
                model.validate_design_row(row);rows.append(row)
    ids=[r['run_fingerprint'] for r in rows]
    if len(set(ids))!=len(ids):raise ValueError('Grid contains duplicate semantic experiments')
    with output.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(model.REQUIRED_COLUMNS)+['fingerprint_version'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    return dict(rows=len(rows),grid=grid,repetitions=repetitions,paired_parameters=sorted(paired),simulation_executed=False)
