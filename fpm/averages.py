"""Reusable graph-first averages; identifiers are loaded as text."""
from pathlib import Path
import pandas as pd
import numpy as np

def summarize(input_path,output,groups,values,graph='graph_sha256'):
    input_path=Path(input_path);output=Path(output)
    if output.exists():raise FileExistsError(output)
    sep='\t' if input_path.suffix.lower()=='.tsv' else ','
    d=pd.read_csv(input_path,sep=sep,dtype=str)
    needed=set(groups+values+[graph])
    if not needed<=set(d):raise ValueError('Missing columns: '+str(sorted(needed-set(d))))
    for value in values:d[value]=pd.to_numeric(d[value],errors='raise')
    if not np.isfinite(d[values].to_numpy()).all():raise ValueError('Filter missing/nonfinite observations explicitly before averaging')
    rows=[]
    for key,cell in d.groupby(groups,dropna=False,sort=True):
        if not isinstance(key,tuple):key=(key,)
        for value in values:
            means=cell.groupby(graph,dropna=False)[value].mean()
            rows.append(dict(zip(groups,key))|dict(observable=value,mean=float(means.mean()),n_graphs=len(means),n_samples=len(cell),graph_mean_sem=float(means.sem()) if len(means)>1 else None,weighting='equal graph weights after within-graph means'))
    output.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(rows).to_csv(output,index=False,float_format='%.17g')
    return dict(rows=len(rows),uncertainty='Graph-mean standard error only; not a confidence interval or calibrated coverage claim')
