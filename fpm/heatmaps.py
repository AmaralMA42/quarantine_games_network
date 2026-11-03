"""Construct measured-cell matrices; missing coordinates remain missing."""
from pathlib import Path
import pandas as pd

def measured_matrix(input_path,output,x,y,value):
    output=Path(output)
    if output.exists():raise FileExistsError(output)
    d=pd.read_csv(input_path,dtype=str)
    for column in [x,y,value]:d[column]=pd.to_numeric(d[column],errors='raise')
    if d.duplicated([x,y]).any():raise ValueError('Filter one topology/observable before constructing a matrix')
    matrix=d.pivot(index=y,columns=x,values=value).sort_index().sort_index(axis=1)
    output.parent.mkdir(parents=True,exist_ok=True);matrix.to_csv(output,float_format='%.17g')
    return dict(measured_cells=len(d),matrix_shape=list(matrix.shape),unmeasured_cells=int(matrix.isna().sum().sum()),interpolation=False)
