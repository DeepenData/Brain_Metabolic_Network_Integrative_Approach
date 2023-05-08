#!/bin/python
"""Genera un grafo a partir de un modelo de COBRA en formato JSON. Lo exporta en
un formato apto para NetworkX y Graph-tool

INPUTS: 
    graphml_path : (str) Direccion del .graphml, eg. "./tmp/graph.graphml"
    nodos_remover : (list) Lista de str de nodos a remover del grafo independientemente. Puede ser lista de listas
OUTPUT: 
    data/centralities.parquet.gzip : (DataFrame) DataFrame de Pandas con centralidades, y multiindice baseline + removidos

"""
# %% --- LOADING MODEL
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
Ggt = gt.load_graph(graphml_path)

nodos_remover = list(Gnx.nodes)[60:70]

# %%
# INICIALIZA LA CONEXION AL SERVIDOR DE RAY
import ray

#if not ray.is_initialized():
#    ray.init()
try:
    ray.init(address='auto', _node_ip_address=os.environ["ip_head"].split(":")[0], _redis_password=os.environ["redis_password"])
except:
    print("No SLURM Cluster detected. Trying a local cluster with default config.")
    ray.init(address='auto',  dashboard_host = '0.0.0.0', _redis_password='5241590000000000')



# Asigns names to the vertexs,
Ggt.vertex_properties["name_ID"] = vertex_names = Ggt.new_vertex_property(
    "string", vals=list(Gnx.nodes)
)


def get_largest_component(grafo):
    import networkx as nx

    largest_component = max(nx.connected_components(grafo), key=len)
    return grafo.subgraph(largest_component)


def compute_centralities_short(graph):
    """Computa las centralidades rapidas y devuelve una DataFrame (no reindexado) con estas"""
    # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
    graph = nx.Graph(graph)
    
    dc = nx.degree_centrality(
        graph
    )  # SANITY CHECK: Deberia ser similar a todo lo demas
    
    #hc = nx.harmonic_centrality(graph)
    ec = nx.eigenvector_centrality(
        graph, max_iter=1000, tol=1e-05, nstart=None, weight=None
    )
    cc = nx.closeness_centrality(graph)
    ic = nx.information_centrality(graph)  # Requiere scipy

    # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
    centralities = {
        "degree_centrality": dc,
        #"harmonic_centrality": hc,
        "eigenvector_centrality": ec,
        "closeness_centrality": cc,
        "information_centrality": ic,
    }

    # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
    centralities = pd.DataFrame(centralities)

    return centralities


def compute_centralities(graph, lite=False, alpha=0.005):
    """Computa las doce centralidades y devuelve una DataFrame (no reindexado) con estas"""
    if lite == False:
        # TODO: supuestamente estos se pueden poner internamente como 'float32', que es suficiente y consume menos memoria
        dc = nx.degree_centrality(
            graph
        )  # SANITY CHECK: Deberia ser similar a todo lo demas
        hc = nx.harmonic_centrality(graph, nbunch=None, distance=None)
        ec = nx.eigenvector_centrality(
            graph, max_iter=1000, tol=1e-05, nstart=None, weight=None
        )
        bc = nx.betweenness_centrality(
            graph, normalized=True, weight=None, endpoints=False, seed=None
        )
        cc = nx.closeness_centrality(graph, distance=None, wf_improved=True)
        lc = nx.load_centrality(graph, cutoff=None, normalized=True, weight=None)
        ic = nx.information_centrality(graph)  # Requiere scipy
        #cbc = nx.communicability_betweenness_centrality(graph)
        kz = nx.katz_centrality(graph, alpha=alpha)
        pr = nx.pagerank(graph, alpha=alpha)
        #acfbc =nx.approximate_current_flow_betweenness_centrality(graph, solver = 'lu',kmax=800)
        #cfbc =nx.current_flow_betweenness_centrality(graph, solver = 'lu')
        #cfcc =nx.current_flow_closeness_centrality(graph, solver = 'lu')
        

        # centralities = [hc, ec, dc, bc, cc, lc, ic, soc, cfcc, cfbc, acfbc, cbc]
        # CREA UN DICCIONARIO DE DICCIONARIOS PARA PASARLE EL OBJETO A PANDAS
        centralities = {
            "degree_centrality": dc,
            "harmonic_centrality": hc,
            "eigenvector_centrality": ec,
            "betweenness_centrality": bc,
            "closeness_centrality": cc,
            "load_centrality": lc,
            "information_centrality": ic,
            #"communicability_betweenness_centrality": cbc,
            "katz_centrality_numpy": kz,
            "pagerank": pr,
        }

        # CONVIERTE LAS CENTRALIDADES A UN DATAFRAME DE PANDAS
        centralities = pd.DataFrame(centralities)

    else:
        # COMPUTA LA VERSION LOW-COST DE LAS CENTRALIDADES
        centralities = compute_centralities_short(graph)

    return centralities


def get_alpha(G):
    """ "Calcula el parametro alpha para PageRank y Katz"""
    tmp = nx.adjacency_matrix(G).toarray()
    tmp = np.linalg.eigvals(tmp)
    tmp = np.max(tmp)
    tmp = np.real(tmp)
    tmp = 0.9 * (1 / tmp)

    return tmp


@ray.remote
def remove_node_get_centrality(graph, removed_node: str, info: bool = False):
    """
    Removes a node from a graph and computes the centralities in it

    [extended_summary]

    Parameters
    ----------
    graph : nx.Graph
        [description]
    node : str
        Key for a node to remove in the graph
    info : bool, optional
        Print extra information, by default False

    Returns
    -------
    removed_node : str
        Key for removed node in the graph
    removed_centrality : pd.DataFrame
        Centralities of every other node removing

    """

    G_removed = graph.copy()  # CREA UNA COPIA DE TRABAJO PARA EL GRAFO

    G_removed.remove_node(removed_node)  # ELIMINA EL NODO A ITERAR
    G_removed = get_largest_component(G_removed)  # ELIMINA NODOS DISCONEXOS

    assert len(graph.nodes) != len(G_removed.nodes), "Ningun nodo removido."
    nodos_removidos = set(graph.nodes) - set(G_removed.nodes)

    if info:
        # INFORMACION EXTRA PARA LOS CALCULOS Y DEMAS
        print(
            f"""
            Nodos originales: {len(graph.nodes)} \n
            Nodos post-remocion: {len(G_removed.nodes)} \n
            Nodos removidos: {len(graph.nodes) - len(G_removed.nodes)} \n
            Removidos: {nodos_removidos} \n
            """
        )

    removed_centrality = compute_centralities(
        G_removed, lite=False, alpha=get_alpha(G_removed)
    )  # CENTRALIDADES

    all_nodes = list(graph.nodes)  # REINDEXANDO PARA INCLUIR REMOVIDO
    removed_centrality = removed_centrality.reindex(all_nodes)

    print(
        f"SANITY_CHECK: {removed_node} \t DIMENSIONES_DATAFRAME: {removed_centrality.shape}"
    )

    removed_centrality.name = str(removed_node)  # RENOMBRADO DEL DATAFRAME

    # SANITY CHECK PARA VER QUE EL DATAFRAME TENGA COMO NAN LOS REMOVIDOS
    for removido in nodos_removidos:
        assert np.isnan(
            removed_centrality.loc[removido, "eigenvector_centrality"]
        ), "Nodo removido tiene un valor no NaN. Error de DataFrame."

    return removed_node, removed_centrality


def remove_several_nodes_get_centrality(graph, nodes_to_remove):
    """Calcula las centralidades para nodos removidos de forma distribuida en un cluster"""

    assert (
        type(nodes_to_remove) == list
    ), "El input de remocion no es una lista. Usa list()"

    # If the node isn't in the model
    short_list = [node for node in nodes_to_remove if node in list(graph.nodes())]

    if len(nodes_to_remove) != len(short_list):
        print(
            f"No todos los nodos estan en la lista ({len(short_list)}/{len(nodes_to_remove)})."
        )

    futures = [remove_node_get_centrality.remote(graph, node) for node in short_list]

    return ray.get(futures)


# %% ---
# DOES THE ACTUAL COMPUTATION STEPS

G = Gnx.copy()

from random import sample
#a sample
#to_remove= sample(nodos_remover,2)
to_remove = nodos_remover

centralidades_perturbadas = remove_several_nodes_get_centrality(G, to_remove)

# TIRA EL CALCULO DE CENTRALIDADES BASE
baseline = compute_centralities(G, lite=False, alpha=get_alpha(G))

# %% ---
# CONCAT PANDAS DATAFRAMES INTO UNIQUE OBJECT WITH MULTIINDEX

final_df = pd.concat(
    [baseline] + [node[1] for node in centralidades_perturbadas],
    axis=0,  # Concate rows
    keys=["baseline"] + [node[0] for node in centralidades_perturbadas],
    names=["removed_node", "metabolite"],
)
# %%
# Saves the dataframe to a parquet file
final_df.to_parquet("data/centralities_long.parquet.gzip")

# %%
