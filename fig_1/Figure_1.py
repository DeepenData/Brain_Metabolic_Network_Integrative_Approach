#%%
import warnings
from graph_tool.topology import extract_largest_component
import pandas as pd
import networkx as nx
import graph_tool as gt
import os
import numpy as np
import networkx as nx
from   cobra.util.array import create_stoichiometric_matrix
import numpy as np
from sklearn.preprocessing import Binarizer
import warnings 
from scipy.sparse import csr_matrix
from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix
import cobra
if "src" in os.getcwd():
    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)
warnings.filterwarnings("ignore")  #
# %%
def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    return grafo.subgraph(largest_component)

def list2attr(grafo, nodos, nombre, atributos):
    """Toma dos listas: nombres de nodos y atributos; y las añade a un grafo

    Parameters
    ----------
    grafo : bipartite_graph
        Un grafo bipartito de NetworkX 
    nodos: list
        Una lista de nombres de nodos
    nombre : str
        Nombre del atributo siendo asignado. Ie. nombre de la columna en Gephi
    atributos : list
        Una lista de valores para los nodos en la lista anterior.
        Ambas listas deben ser del mismo largo. 
    
    Returns
    -------
    grafo: bipartite_graph
        Un grafo bipartito con un nuevo atributo para un set de nodos "nodos". 
    """
    assert len(nodos) == len(atributos), "Ambas listas deben ser del mismo largo."
    tmp_list = { nodos[i] : atributos[i] for i in range(0, len(atributos)) }
    
    from networkx import set_node_attributes
    set_node_attributes(grafo, tmp_list, nombre ) # Añade los atributos
    
    return grafo
# %%
#
model    = cobra.io.load_json_model('data/models/stimulated_NAMU.json')
S_matrix = create_stoichiometric_matrix(model)
S_matrix = abs(S_matrix)
S_matrix = Binarizer().fit_transform(S_matrix).astype(int)
grafo    = from_biadjacency_matrix(csr_matrix(S_matrix)) # Usa la dispersa para un grafo bipartito
#
metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
metabolites_n     = len(model.metabolites) # Numero de metabolitos
metabolites_names = [model.metabolites[i].id for i in range(metabolites_n)]
#
reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones
reactions_n       = len(model.reactions)   # Numero de reacciones
reactions_names   = [model.reactions[i].id for i in range(reactions_n)]
#
names_mapped      =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reactions_names))
grafo             =  nx.relabel_nodes(grafo, names_mapped)
# %%

def pseudoLog10(x):
    return np.arcsinh(x/2) / np.log(10)

def scale_data(x):   

    from sklearn.preprocessing import MinMaxScaler
    scaler       = MinMaxScaler(feature_range=(0, 15))  
    return scaler.fit_transform( pseudoLog10(abs(x)).reshape(-1, 1) ).astype(np.float32).flatten()     #queda como un array (pierde el multiindice)

FBA_step01_FBA_solution = pd.read_parquet("data/FBA_step01_FBA_solution.parquet.gzip")
with_cells              = pd.read_parquet("data/FBA_step01_FBA_solution_with_cells.parquet.gzip")

  
flux_arr = np.array(FBA_step01_FBA_solution.fluxes.values).reshape(-1, 1)
fluxes   =  scale_data(flux_arr)
#fluxes = abs(np.log(abs(np.array(FBA_step01_FBA_solution.fluxes.values).reshape(-1, 1))))
#fluxes = np.nan_to_num(fluxes, copy=True, posinf=0, neginf=0)


grafo                   = list2attr(grafo, FBA_step01_FBA_solution.index.values, 'Flux', 
                                    list(fluxes))



sensitivity_arr = np.array(FBA_step01_FBA_solution.reduced_costs.values).reshape(-1, 1)
sensitivity     = scale_data(sensitivity_arr)
grafo           = list2attr(grafo, FBA_step01_FBA_solution.index.values, 'Sensitivity', 
                                    list(sensitivity) )


cell = list(pd.Series(with_cells.cell.values).replace('Astrocyte', 0, regex=False).replace('Neuron',    1, regex=False))



grafo                   = list2attr(grafo,with_cells.ID.values , 'Cell', cell)
G = grafo.copy()
G.remove_nodes_from(pd.Series(reactions_names).str.extractall('^(DM_.*|.*sink.*|EX_.*)').loc[:,0].values)
G = get_largest_component(G)
print(
len(list(grafo.nodes)),
len(list(G.nodes)))

shapes = \
dict(zip(list(G.nodes),
1*(np.array(list(nx.get_node_attributes(G, "bipartite").values()))==0) ))


cells = nx.get_node_attributes(G, "Cell")
cells.update(dict(zip(metabolites_names, [2]*len(metabolites_names))))


nx.set_node_attributes(G, shapes, "shapes")
nx.set_node_attributes(G, cells, "cells")


nx.write_graphml(G, "data/models/pruned_stimulated_w_attrs.graphml")


# %%
from graph_tool import *
from graph_tool.draw import *
import matplotlib

GT    = load_graph("data/models/pruned_stimulated_w_attrs.graphml")

np.random.seed(23)
seed_rng(42)
pos = sfdp_layout(GT,groups=GT.vertex_properties['cells'],max_iter =0,mu =2.6, gamma = 1.2, C = 4.15, p = 3)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap



colors =["green", "midnightblue", 
         "maroon"]

newcmp= matplotlib.colors.ListedColormap(colors)

graph_draw(GT, pos=pos,vertex_font_size = 7,vcmap=newcmp, output_size=  (450,400),  \
  vertex_fill_color = GT.vertex_properties['cells'],
    #vertex_halo = True,vertex_halo_size = .8, #vertex_halo_color = 'purple' ,
   output= "fig_1/pruned_stimulated_w_attrs.svg",
   vertex_shape = GT.vertex_properties['shapes'],
           vertex_size= 3.9, bg_color = None, edge_pen_width=0.4, edge_color = "grey",
            edge_marker_size = 3, fmt = 'svg')
# %%
graph_draw(GT,pos=pos,vertex_font_size = 7,vcmap=newcmp ,output_size= (700,700),  \
   vertex_size =  GT.vertex_properties['Flux'],
    vertex_fill_color = GT.vertex_properties['cells'],
        vertex_halo = True,vertex_halo_size = 1.3, vertex_halo_color = 'darkslategray' ,
   #vertex_text =  GT.vertex_properties['_graphml_vertex_id'],\
   output= "fig_1/pruned_stimulated_w_attrs_Flux.svg",
   vertex_shape = GT.vertex_properties['shapes'],
            bg_color = None, edge_pen_width=0.3, edge_color = "grey",
            edge_marker_size = 0, fmt = 'svg')
# %%
graph_draw(GT,pos=pos,vertex_font_size = 7,vcmap=newcmp ,output_size= (700,700),  \
   vertex_size =  GT.vertex_properties['Sensitivity'],
    vertex_fill_color = GT.vertex_properties['cells'],
    vertex_halo = True,vertex_halo_size = 1.3, vertex_halo_color = 'darkslategray' ,
   #vertex_text =  GT.vertex_properties['_graphml_vertex_id'],\
   output= "fig_1/pruned_stimulated_w_attrs_Sensitivity.svg",
   vertex_shape = GT.vertex_properties['shapes'],
            bg_color = None, edge_pen_width=.3, edge_color = "grey",
            edge_marker_size = 0, vertex_anchor=0.8, fmt = 'svg')

# %%