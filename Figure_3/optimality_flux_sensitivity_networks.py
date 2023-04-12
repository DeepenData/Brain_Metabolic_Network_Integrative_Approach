# %% 
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
figure1_path = local_working_directory + "/Figure_1"
figure2_path = local_working_directory + "/Figure_2"
figure3_path = local_working_directory + "/Figure_3"
from aa_functions import list2attr
import warnings 
import pandas as pd
import networkx as nx

#reaction_projected_optimality_weighted
reaction_projected_stimnulated = nx.read_gpickle( figure2_path + "/reaction_projected_stimnulated.gpickle")
optimality_values = pd.read_csv(figure3_path + '/optimality_values.csv')
reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = list(optimality_values.ID ), \
                                  nombre = "optimality", atributos = optimality_values.total)

nx.write_gexf(reaction_projected_stimnulated,  figure3_path + "/reaction_projected_flux_sensi_centr_opti_weighted.gexf")
#reaction_projected_stimnulated.nodes(data=True)    

nx.write_gpickle(reaction_projected_stimnulated,  figure3_path + "/reaction_projected_flux_sensi_centr_opti_weighted.gpickle")

# %% 

stimulated_bipartite_graph = nx.read_gpickle( figure1_path + "/stimulated_bipartite_graph.gpickle")

stimulated_bipartite_graph  = list2attr(stimulated_bipartite_graph, nodos = list(optimality_values.ID ), \
                                  nombre = "Optimality", atributos = optimality_values.total)

stimulated_bipartite_graph  = list2attr(stimulated_bipartite_graph, nodos = list(optimality_values.ID ), \
                                  nombre = "Flux", atributos = optimality_values.Flux)

stimulated_bipartite_graph  = list2attr(stimulated_bipartite_graph, nodos = list(optimality_values.ID ), \
                                  nombre = "Sensitivity", atributos = optimality_values.Sensitivity)                                                                    

nx.write_gexf(stimulated_bipartite_graph,  figure3_path + "/bipartite_flux_sensi_opti_weighted.gexf")
#   stimulated_bipartite_graph.nodes(data=True) 
# %%                               