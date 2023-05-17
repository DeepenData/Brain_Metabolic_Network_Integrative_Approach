# %%
import ray

import networkx as nx
import numpy    as np
import pandas   as pd 
import inspect
import warnings 
import copy
warnings.filterwarnings("ignore")

def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = copy.deepcopy(grafo.subgraph(largest_component))
    assert nx.is_connected(G)
    return G


@ray.remote # TODO: incluir parametros del metodo
def compute_centralities(graph, alpha=0.005):
    """Computa las doce centralidades y devuelve una DataFrame (no reindexado) con estas"""
    #if lite == False:
    # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
    dc    = nx.degree_centrality(graph) # SANITY CHECK: Deberia ser similar a todo lo demas
    hc    = nx.harmonic_centrality(graph, nbunch=None, distance=None)
    ec    = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05, nstart=None, weight=None)
    bc    = nx.betweenness_centrality(graph, normalized=True, weight=None, endpoints=False, seed=None)
    cc    = nx.closeness_centrality(graph, distance=None, wf_improved=True)
    lc    = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
    ic    = nx.information_centrality(graph) # Requiere scipy
    cbc   = nx.communicability_betweenness_centrality(graph) 
    kz    = nx.katz_centrality(graph, alpha = alpha)
    pr    = nx.pagerank(graph, alpha = alpha)
    
    cfcc  = nx.current_flow_closeness_centrality(graph)
    cfbc  = nx.current_flow_betweenness_centrality(graph)
    #soc   = nx.second_order_centrality(graph)
    lac   = nx.laplacian_centrality(graph)
    

    #centralities = [hc, ec, dc, bc, cc, lc, ic, soc, cfcc, cfbc, acfbc, cbc]
    # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
    centralities = {
        'degree_centrality' : dc ,
        'harmonic_centrality' : hc ,
        'eigenvector_centrality' : ec ,
        'betweenness_centrality' : bc ,
        'closeness_centrality' : cc ,
        'load_centrality' : lc ,
        'information_centrality' : ic ,
        'communicability_betweenness_centrality' : cbc ,
        'katz_centrality_numpy' : kz ,
        'pagerank' : pr,
        'current_flow_closeness_centrality':cfcc,
        'current_flow_betweenness_centrality':cfbc,
        'laplacian_centrality': lac
        }
    
    # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
    centralities = pd.DataFrame( centralities )
    # else: 
    #     # COMPUTA LA VERSION LOW-COST DE LAS CENTRALIDADES  
    #     #centralities = compute_centralities_short(graph)
    return centralities
# %%
G          = nx.Graph(nx.read_graphml("tests/secretome/G.graphml"))
names_list = [data['name'] for _, data in G.nodes(data=True)]



mapping = dict(zip(G.nodes, names_list))

G = nx.relabel_nodes(G, mapping, copy=False)
G.nodes(data=True)
# %%


G2 = get_largest_component(G)


g2_df = ray.get( compute_centralities.remote(G2) )
g2_df.to_parquet("tests/secretome/centralities.parquet")


import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import matplotlib.pyplot as plt

df = pd.read_parquet("tests/secretome/centralities.parquet")

# Standardizing rows to have zero mean and unit variance
standard_scaler = StandardScaler()
df_standard_scaled = pd.DataFrame(standard_scaler.fit_transform(df), columns=df.columns, index=df.index)

# Scaling standardized rows to fall into the range [0, 1]
min_max_scaler = MinMaxScaler()
df_scaled = pd.DataFrame(min_max_scaler.fit_transform(df_standard_scaled), columns=df.columns, index=df.index)

# Creating the clustermap
def draw_clustermap(df, colorbar_pos):
    # Creating the clustermap
    g = sns.clustermap(df,  metric='euclidean', method='centroid')

    # Adjust colorbar position
    g.cax.set_position(colorbar_pos)

    # Add label to the colorbar
    g.cax.yaxis.set_label_position("right")
    g.cax.set_ylabel('Centrality', rotation=90, va="bottom", labelpad=15)
    
        # Rotate xticklabels and adjust size
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=10)

    # Add padding
    g.ax_heatmap.tick_params(axis='x', pad=2)


# Example usage
draw_clustermap(df_scaled, [1., 0.36, 0.03, 0.44])  # [left, bottom, width, height]


#sns.clustermap(df_scaled, metric='euclidean', method='centroid')
