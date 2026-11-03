"""Reconstitute original M25 edges using recorded seeds and original serialization."""
from pathlib import Path
import hashlib,json
import pandas as pd,networkx as nx
R=Path(__file__).resolve().parent
assert nx.__version__=='3.6.1',nx.__version__
data=pd.read_csv(R/'inputs/m25_metadata/M25_RRG_GRAPH_BANK_AUDIT.csv')
data=data[data.log2N_nominal==14];assert len(data)==12
O=R/'results/reconstructed_fig10_graphs';O.mkdir(parents=True,exist_ok=True)
rows=[]
for t in data.itertuples():
    g=nx.random_regular_graph(int(t.degree),int(t.population),seed=int(t.seed32))
    # The historical file follows NetworkX native node/adjacency iteration.
    # Sorting the edges changes file bytes, despite an identical contact graph.
    edges=list(g.edges());b=(f'{len(g)} {len(edges)}\n'+''.join(f'{u} {v}\n' for u,v in edges)).encode('ascii')
    h=hashlib.sha256(b).hexdigest();assert h==t.sha256,(t.topology,t.graph_index,h,t.sha256)
    name=f'{t.topology}_N2p14_g{t.graph_index:02d}.edges';(O/name).write_bytes(b)
    rows.append(dict(file=name,sha256=h,seed32=int(t.seed32),byte_identical_to_original=True))
(R/'qa/FIG10_GENERATOR_REPLAY.json').write_text(json.dumps(dict(status='PASS',generator='networkx.random_regular_graph',version=nx.__version__,original_graphs_reconstituted=12,new_graph_samples=0,serialization='native nx.Graph.edges iteration; LF ASCII',files=rows),indent=2)+'\n',encoding='utf-8')
print('Original graph generator replay:12/12 byte-identical')
