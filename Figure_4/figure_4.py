# %% 
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

reaction_projected_stimnulated = nx.read_gpickle( figure3_path + "/reaction_projected_flux_sensi_centr_opti_weighted.gpickle")
hubs_info = pd.read_csv(figure4_path + '/hubs_info.csv')
node_list0 = list(reaction_projected_stimnulated.nodes()) 
node_list = list(map(lambda st: str.replace(st, "_", ""), node_list0))
new_node_list = list(map(lambda st: str.replace(st, ",", ""), node_list))
mapping = dict(zip(node_list0, new_node_list))
reaction_projected_stimnulated= nx.relabel_nodes(reaction_projected_stimnulated, mapping) 
# %% 
reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = list(hubs_info.ID ), \
                                  nombre = "module", atributos = hubs_info.Module)

nx.write_gexf(reaction_projected_stimnulated,  figure4_path + "/reaction_projected_flux_sensi_centr_opti_weighted.gexf")
#reaction_projected_stimnulated.nodes(data=True)    

nx.write_gpickle(reaction_projected_stimnulated,  figure4_path + "/reaction_projected_flux_sensi_centr_opti_weighted.gpickle")
