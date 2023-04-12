# %% 
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
figure1_path = local_working_directory + "/Figure_1"
figure2_path = local_working_directory + "/Figure_2"
figure3_path = local_working_directory + "/Figure_3"
figure4_path = local_working_directory + "/Figure_4"
supp_path    = local_working_directory + "/Supplementary/"
supp_path_results    = local_working_directory + "/Supplementary/results/"
from aa_functions import *
from networkx import bipartite
from graph_tool import *
#Random
#Creación de dos grafos random (generate_rbd_graph): G1 y H1. 
one_core_diRandom, two_core_diRandom = get_two_core_directed_random(A = 20, B = 16, p = 0.10)
#Free-scale
#Generar una secuencia de grados (A_degree_sequence) en base al grafo anterior G1 para crear una grafo con "Prefential Attachment" 
#que da origen a la topología de 'escala-libre'.
#crear A_degree_sequence:
degree_sequence = [d for n, d in one_core_diRandom.degree()]
A, B = bipartite.sets(one_core_diRandom)
A_degree_sequence = degree_sequence[0:len(A)]
one_core_diFreeScale, two_core_diFreeScale = \
    get_two_core_directed_FreeScale(A = A_degree_sequence, p = 0.05)
#Projections
Proj_one_core_diRandom = make_rxn_projection(one_core_diRandom)
Proj_two_core_diRandom = make_rxn_projection(two_core_diRandom)
Proj_one_core_diFreeScale= make_rxn_projection(one_core_diFreeScale)
Proj_two_core_diFreeScale= make_rxn_projection(two_core_diFreeScale)
import networkx as nx
nx.write_graphml(one_core_diRandom, supp_path_results + "one_core_diRandom.graphml")
nx.write_graphml(Proj_one_core_diRandom, supp_path_results + "Proj_one_core_diRandom.graphml")
nx.write_graphml(two_core_diRandom, supp_path_results + "two_core_diRandom.graphml")
nx.write_graphml(Proj_two_core_diRandom, supp_path_results + "Proj_two_core_diRandom.graphml")
nx.write_graphml(one_core_diFreeScale, supp_path_results + "one_core_diFreeScale.graphml")
nx.write_graphml(Proj_one_core_diFreeScale, supp_path_results + "Proj_one_core_diFreeScale.graphml")
nx.write_graphml(two_core_diFreeScale, supp_path_results + "two_core_diFreeScale.graphml")
nx.write_graphml(Proj_two_core_diFreeScale, supp_path_results + "Proj_two_core_diFreeScale.graphml")

set_subgraphs_attr_bipartite( supp_path_results + "one_core_diRandom.graphml", supp_path_results +'one_core_diRandom.png')
set_subgraphs_attr( supp_path_results + "Proj_one_core_diRandom.graphml", supp_path_results +'Proj_one_core_diRandom.png')
set_subgraphs_attr_bipartite( supp_path_results + "two_core_diRandom.graphml", supp_path_results +'two_core_diRandom.png')
set_subgraphs_attr( supp_path_results + "Proj_two_core_diRandom.graphml", supp_path_results +'Proj_two_core_diRandom.png')
set_subgraphs_attr_bipartite( supp_path_results + "one_core_diFreeScale.graphml", supp_path_results +'one_core_diFreeScale.png')
set_subgraphs_attr( supp_path_results + "Proj_one_core_diFreeScale.graphml", supp_path_results +'Proj_one_core_diFreeScale.png')
set_subgraphs_attr_bipartite( supp_path_results + "two_core_diFreeScale.graphml", supp_path_results +'two_core_diFreeScale.png')
set_subgraphs_attr( supp_path_results + "Proj_two_core_diFreeScale.graphml", supp_path_results +'Proj_two_core_diFreeScale.png')
# %%
