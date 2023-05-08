# %%
import warnings
from graph_tool.topology import extract_largest_component
from networkx.algorithms.centrality import current_flow_betweenness

warnings.filterwarnings("ignore")  #
import pandas as pd
import networkx as nx
import graph_tool.centrality as gt_cent
import graph_tool.clustering as gt_clus
import graph_tool as gt
import os
import numpy as np


if "src" in os.getcwd():
    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)
    
graphml_path: str = "./tmp/graph.graphml"
Gnx = nx.read_graphml(graphml_path)  #
Gnx = nx.Graph(Gnx)

# %%
G = Gnx.copy()
n = len(G.nodes)
m = G.size()

average_degree = 2*m/n
density        = m/(n*(n-1)/2)

print("number of edges" , m)
print("number of nodes" ,n)
print("average_degree" , average_degree)
print("density" ,density)
# %%
#Distance based measures 


