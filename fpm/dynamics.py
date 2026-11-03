"""One execution interface to the retained, generic Quarantine Games kernel.

Epidemic transition rates remain in engine/src/fpm_network_sim.cpp. This module
builds that source, checks run/graph identities and runs parameterized designs.
"""
from pathlib import Path
import csv,hashlib,json,os,subprocess
from .model import validate_design_row
ROOT=Path(__file__).resolve().parents[1]

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def build_engine(compiler):
    target=ROOT/'build'/('fpm_network_sim.exe' if os.name=='nt' else 'fpm_network_sim')
    target.parent.mkdir(exist_ok=True)
    command=[str(compiler),'-std=c++17','-O2','-I'+str(ROOT/'engine/include'),str(ROOT/'engine/src/fpm_network_sim.cpp'),'-o',str(target)]
    if os.name=='nt':command[1:1]=['-static','-static-libgcc','-static-libstdc++']
    result=subprocess.run(command,capture_output=True,text=True,timeout=300)
    if result.returncode:raise RuntimeError(result.stderr[-4000:])
    sources={p.relative_to(ROOT).as_posix():digest(p) for p in (ROOT/'engine').rglob('*') if p.is_file()}
    receipt=dict(binary=target.name,binary_sha256=digest(target),source_files=sources,command=command)
    (ROOT/'build/engine_build.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    return receipt

def read_design(path):
    path=Path(path).resolve()
    with path.open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows:raise ValueError('Empty design')
    for row in rows:
        validate_design_row(row)
        if row['topology']=='well_mixed':continue
        graph=Path(row['graph_path'])
        if not graph.is_absolute():graph=path.parent/graph
        if digest(graph)!=row['graph_sha256']:raise ValueError('Graph checksum mismatch')
        if row['graph_certificate_path']!='-':
            cert=Path(row['graph_certificate_path'])
            if not cert.is_absolute():cert=path.parent/cert
            if digest(cert)!=row['graph_certificate_sha256']:raise ValueError('Periodic certificate checksum mismatch')
    ids=[x['run_fingerprint'] for x in rows]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate planned identities')
    return rows

def run_design(design,output,backend='fenwick_exact'):
    design=Path(design).resolve();output=Path(output).resolve()
    if output.exists():raise FileExistsError('Choose a new output path; existing results are preserved')
    rows=read_design(design)
    receipt=json.loads((ROOT/'build/engine_build.json').read_text())
    binary=ROOT/'build'/receipt['binary']
    if digest(binary)!=receipt['binary_sha256']:raise ValueError('Engine binary changed after build')
    for path,sha in receipt['source_files'].items():
        if digest(ROOT/path)!=sha:raise ValueError('Rebuild after changing the shared engine source')
    output.parent.mkdir(parents=True,exist_ok=True)
    result=subprocess.run([str(binary),str(design),str(output),backend],cwd=design.parent,capture_output=True,text=True,timeout=3600)
    if result.returncode:raise RuntimeError(result.stderr[-4000:])
    with output.open(newline='',encoding='utf-8') as f:observed=list(csv.DictReader(f,delimiter='\t'))
    if [x['run_fingerprint'] for x in observed]!=[x['run_fingerprint'] for x in rows]:raise ValueError('Planned and produced identities differ')
    record=dict(kernel_sha256=digest(binary),design_sha256=digest(design),output_sha256=digest(output),runs=len(rows),backend=backend)
    output.with_suffix(output.suffix+'.provenance.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    return record
