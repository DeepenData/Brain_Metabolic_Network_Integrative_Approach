# %% 

import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
figure_path = local_working_directory + "/Figure_2"
from aa_functions import cobra_to_bipartite, get_largest_component, list2attr
import warnings 
import pandas as pd
import networkx as nx

warnings.filterwarnings("ignore")
nodes_centrality = pd.read_csv(figure_path + "/total_centrality_by_node.csv")

reaction_projected_stimnulated = nx.read_gpickle(local_working_directory + "/Figure_1/"  + "reaction_projected_stimnulated.gpickle")
reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = nodes_centrality.ID , \
                                            nombre = "Centrality",   atributos = nodes_centrality.total_centrality)


#reaction_projected_stimnulated.nodes(data = True)
nx.write_gpickle(reaction_projected_stimnulated, figure_path + "/reaction_projected_stimnulated.gpickle")

nx.write_gexf(reaction_projected_stimnulated,  figure_path + "/reaction_projected_stimnulated.gexf")