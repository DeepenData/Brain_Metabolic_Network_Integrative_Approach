# %% Importa el modelo JSON
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
figure1_path = local_working_directory + "/Figure_1"
figure2_path = local_working_directory + "/Figure_2"
figure3_path = local_working_directory + "/Figure_3"
figure4_path = local_working_directory + "/Figure_4"
from aa_functions import list2attr
import warnings 
import pandas as pd
import networkx as nx

reaction_projected = nx.read_gpickle(figure3_path + "/reaction_projected_flux_sensi_centr_opti_weighted.gpickle")
optimality_values  = pd.read_csv(figure3_path + '/optimality_values.csv')


reaction_projected  = list2attr(reaction_projected, nodos = list(optimality_values.ID ), \
                                  nombre = "Optimality", atributos = optimality_values.total)

reaction_projected  = list2attr(reaction_projected, nodos = list(optimality_values.ID ), \
                                  nombre = "Flux", atributos = optimality_values.Flux)

reaction_projected  = list2attr(reaction_projected, nodos = list(optimality_values.ID ), \
                                  nombre = "Sensitivity", atributos = optimality_values.Sensitivity)                                                                    



Astrocyte_nodes = [n for n, d in reaction_projected.nodes(data=True) if d["type"] == "Astrocyte"] 
Neuron_nodes    = [n for n, d in reaction_projected.nodes(data=True) if d["type"] == "Neuron"] 


Neuron    = reaction_projected.subgraph(Neuron_nodes)
Astrocyte = reaction_projected.subgraph(Astrocyte_nodes)
# %%
nx.write_gexf(Neuron,  figure4_path + "/Neuron.gexf")
nx.write_gexf(Astrocyte,  figure4_path + "/Astrocyte.gexf")


nx.write_gpickle(Neuron,  figure4_path + "/Neuron.gpickle")
nx.write_gpickle(Astrocyte,  figure4_path + "/Astrocyte.gpickle")