import warnings
warnings.filterwarnings("ignore")  #
import logging
import numpy as np
import networkx as nx
import pandas as pd
import ray

# Define the format
log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

# Configure the basic logger
logging.basicConfig(filename='centralities.log', level=logging.INFO, format=log_format, datefmt='%m/%d/%Y %I:%M:%S %p')


def get_largest_connected_component(graph: nx.Graph) -> nx.Graph:
    """
    Find and return the largest connected component of a graph.

    Args:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The largest connected subgraph.
    """
    assert isinstance(graph, nx.Graph), "Input must be a networkx Graph."
    
    try:
        largest_component = max(nx.connected_components(graph), key=len)
        largest_connected_subgraph = graph.subgraph(largest_component).copy()
        assert nx.is_connected(largest_connected_subgraph), "Output subgraph is not connected."
    except Exception as e:
        logging.error(f"Failed to compute largest connected component: {str(e)}")
        raise

    return largest_connected_subgraph



def calculate_quick_centralities(graph: nx.Graph) -> pd.DataFrame:
    """
    Compute some of the faster-to-calculate centralities and return them in a DataFrame.

    Args:
        graph (nx.Graph): The input graph.

    Returns:
        pd.DataFrame: DataFrame with the calculated centralities.
    """
    assert isinstance(graph, nx.Graph), "Input must be a networkx Graph."

    try:
        # Degree centrality
        degree_centrality = nx.degree_centrality(graph)

        # Eigenvector centrality
        eigenvector_centrality = nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05)

        # Closeness centrality
        closeness_centrality = nx.closeness_centrality(graph)

        # Information centrality
        info_centrality = nx.information_centrality(graph)

        centralities = {
            "degree_centrality": degree_centrality,
            "eigenvector_centrality": eigenvector_centrality,
            "closeness_centrality": closeness_centrality,
            "information_centrality": info_centrality,
        }

        df_centralities = pd.DataFrame(centralities)
    except Exception as e:
        logging.error(f"Failed to compute quick centralities: {str(e)}")
        raise

    return df_centralities

def calculate_centralities(graph: nx.Graph, use_quick_measurements: bool=False, alpha: float=0.005) -> pd.DataFrame:
    """
    Compute several centralities and return them in a DataFrame.

    Args:
        graph (nx.Graph): The input graph.
        use_quick_measurements (bool): Flag to use quick measurements. Default is False.
        alpha (float): The damping factor for Katz centrality and PageRank. Default is 0.005.

    Returns:
        pd.DataFrame: DataFrame with the calculated centralities.
    """
    assert isinstance(graph, nx.Graph), "Input must be a networkx Graph."
    assert isinstance(use_quick_measurements, bool), "use_quick_measurements must be a boolean."
    assert isinstance(alpha, (int, float)), "alpha must be a numeric type."

    try:
        if not use_quick_measurements:
            centralities = {
                "degree_centrality": nx.degree_centrality(graph),
                "harmonic_centrality": nx.harmonic_centrality(graph),
                "eigenvector_centrality": nx.eigenvector_centrality(graph, max_iter=1000, tol=1e-05),
                "betweenness_centrality": nx.betweenness_centrality(graph, normalized=True),
                "closeness_centrality": nx.closeness_centrality(graph, wf_improved=True),
                "load_centrality": nx.load_centrality(graph, normalized=True),
                "information_centrality": nx.information_centrality(graph),
                "katz_centrality": nx.katz_centrality(graph, alpha=alpha),
                "pagerank": nx.pagerank(graph, alpha=alpha),
            }
        else:
            centralities = calculate_quick_centralities(graph)
    except Exception as e:
        logging.error(f"Failed to compute centralities: {str(e)}")
        raise

    return pd.DataFrame(centralities)


def compute_alpha_for_graph(graph: nx.Graph) -> float:
    """
    Calculate the alpha parameter for Katz centrality and PageRank based on the graph's largest eigenvalue.

    Args:
        graph (nx.Graph): The input graph.

    Returns:
        float: The calculated alpha.
    """
    assert isinstance(graph, nx.Graph), "Input must be a networkx Graph."

    try:
        adjacency_matrix = nx.adjacency_matrix(graph).toarray()
        largest_eigenvalue = np.max(np.linalg.eigvals(adjacency_matrix))
        alpha = 0.9 * (1 / np.real(largest_eigenvalue))
    except Exception as e:
        logging.error(f"Failed to compute alpha: {str(e)}")
        raise

    return alpha



@ray.remote
def remove_node_and_calculate_centralities(graph: nx.Graph, node_to_remove, verbose: bool=False):
    """
    Remove a node from the graph, then calculate and return the centralities of the remaining nodes.

    Args:
        graph (nx.Graph): The input graph.
        node_to_remove: The node to be removed.
        verbose (bool): Flag to print verbose messages. Default is False.

    Returns:
        tuple: A tuple containing the removed node and a DataFrame of the new centralities.
    """
    assert isinstance(graph, nx.Graph), "Input must be a networkx Graph."
    assert isinstance(verbose, bool), "verbose must be a boolean."

    try:
        graph_copy = graph.copy()
        graph_copy.remove_node(node_to_remove)
        largest_connected_subgraph = get_largest_connected_component(graph_copy)
    except Exception as e:
        logging.error(f"Failed to compute largest connected component: {str(e)}")
        raise

    assert len(graph.nodes) != len(largest_connected_subgraph.nodes), "No node was removed."

    removed_nodes = set(graph.nodes) - set(largest_connected_subgraph.nodes)

    if verbose:
        logging.info(f"Nodes removed: {removed_nodes}")

    try:
        new_centralities = calculate_centralities(largest_connected_subgraph, use_quick_measurements=False, alpha=compute_alpha_for_graph(largest_connected_subgraph))
        new_centralities = new_centralities.reindex(list(graph.nodes))
        new_centralities.name = str(node_to_remove)
    except Exception as e:
        logging.error(f"Failed to calculate centralities: {str(e)}")
        raise

    # Check that removed nodes have NaN centralities
    for removed_node in removed_nodes:
        assert np.isnan(new_centralities.loc[removed_node, "eigenvector_centrality"]), "Removed node has non-NaN value."

    return node_to_remove, new_centralities



def remove_nodes_and_calculate_centralities(graph: nx.Graph, nodes_to_remove: list):
    """
    Remove a list of nodes from the graph and calculate the centralities after each removal.

    Args:
        graph (nx.Graph): The input graph.
        nodes_to_remove (list): List of nodes to be removed.

    Returns:
        list: List of tuples with removed node and new centralities.
    """
    # Checking input types
    assert isinstance(graph, nx.Graph), "Input must be a networkx Graph."
    assert isinstance(nodes_to_remove, list), "nodes_to_remove must be a list."

    # Check if all nodes to remove are in the graph, and create a list of nodes in the graph
    nodes_in_graph = [node for node in nodes_to_remove if node in graph.nodes()]

    # If any nodes_to_remove were not found in the graph, log a warning
    if len(nodes_to_remove) != len(nodes_in_graph):
        logging.warning(f"Not all nodes are in the graph ({len(nodes_in_graph)}/{len(nodes_to_remove)}).")

    try:
        # Here, we use Ray to execute the remove_node_and_calculate_centralities function in parallel for each node
        # The .remote() function call tells Ray to execute the function as a remote task
        # This line does not actually execute the tasks yet, but creates futures for them
        futures = [remove_node_and_calculate_centralities.remote(graph, node) for node in nodes_in_graph]

        # This line triggers the execution of the remote tasks and blocks until all tasks have completed
        # The results of the tasks are returned as a list in the same order as the futures
        result = ray.get(futures)
    except Exception as e:
        # If anything goes wrong during the execution of the tasks, log an error and re-raise the exception
        logging.error(f"Failed to remove nodes and calculate centralities: {str(e)}")
        raise

    # Return the result, which is a list of (removed_node, new_centralities) tuples
    return result


if __name__ == '__main__':
    graphml_path: str = "graph.graphml"
    Gnx = nx.read_graphml(graphml_path)  #
    Gnx = nx.Graph(Gnx)

    G = Gnx.copy()

    to_remove = list(Gnx.nodes)[68:70]


    ray.init(ignore_reinit_error=True)


    centralidades_perturbadas = remove_nodes_and_calculate_centralities(G, to_remove)
    baseline                  = calculate_centralities(G)

    ray.shutdown()
    
    import pickle
    import os

    filename = 'centralities.pickle'

    # Make sure the necessary variables are defined
    assert 'centralidades_perturbadas' in globals(), "centralidades_perturbadas not defined"
    assert 'baseline' in globals(), "baseline not defined"

    try:
        with open(filename, 'wb') as handle:
            pickle.dump({
                'centralidades_perturbadas': centralidades_perturbadas,
                'baseline': baseline
            }, handle, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as e:
        logging.error(f"Error while saving the dictionary: {e}")
    else:
        if os.path.exists(filename) and os.path.getsize(filename) > 0:
            logging.info("The dictionary was successfully saved to 'centralities.pickle'.")
        else:
            logging.warning("The dictionary was not saved correctly.")

    try:
        logging.info(f'''
            total nodes: {centralidades_perturbadas[0][1].shape[0]}
            total centralites: {baseline.shape[1]}
            nodes removed: {len(centralidades_perturbadas)}
            ''')
    except Exception as e:
        logging.error(f"Error while logging: {e}")
    # Create a file handler
    
    # Set level of logging

# with open('centralities.pickle', 'rb') as handle:
#     data = pickle.load(handle)

#     centralidades_perturbadas = data['centralidades_perturbadas']
#     baseline = data['baseline']


#ray up       aws-ray-cluster.yml
#ray rsync-up aws-ray-cluster.yml graph.graphml graph.graphml                
#ray submit   aws-ray-cluster.yml centralities.py



