# %% 
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
figure_path = local_working_directory + "/Figure_1"
from aa_functions import cobra_to_bipartite, get_largest_component, list2attr,solution2attr, cobra_to_networkx_rxn_projection
import cobra
import warnings 
import pandas as pd
import networkx as nx
warnings.filterwarnings("ignore")
stimulated      = cobra.io.load_json_model(local_working_directory + '/stimulated_2021.json')
nodes_attr_df   = pd.read_csv(local_working_directory + '/node_attributes.csv')
stimulated_graph = cobra_to_bipartite(stimulated)

fba_solution        = stimulated.optimize()
stimulated_graph    = solution2attr(fba_solution, stimulated_graph, estandarizar=True)

rxns   =  [n for n, d in stimulated_graph.nodes(data=True) if d["bipartite"] == 1] 
stimulated_graph  = list2attr(stimulated_graph, nodos = nodes_attr_df.ID , nombre = "type", atributos = nodes_attr_df.Node)
stimulated_graph = get_largest_component(stimulated_graph)
nx.write_gexf(stimulated_graph,figure_path+ "/stimulated_bipartite_graph.gexf")
nx.write_gpickle(stimulated_graph,figure_path+ "/stimulated_bipartite_graph.gpickle")
# %%
reaction_projected_stimnulated = cobra_to_networkx_rxn_projection(stimulated)
fluxes        = fba_solution.fluxes.to_frame()

reduced_costs = fba_solution.reduced_costs.to_frame()
reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = fluxes.index , \
                                  nombre = "fluxes", atributos = fluxes.fluxes)
reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = reduced_costs.index , \
                                  nombre = "reduced_costs", atributos = reduced_costs.reduced_costs)
reaction_projected_stimnulated  = list2attr(reaction_projected_stimnulated, nodos = nodes_attr_df.ID , \
    nombre = "type", atributos = nodes_attr_df.Node)
reaction_projected_stimnulated = get_largest_component(reaction_projected_stimnulated)

nx.write_gpickle(reaction_projected_stimnulated, figure_path + "/reaction_projected_stimnulated.gpickle")
nx.write_gexf(reaction_projected_stimnulated,  figure_path + "/reaction_projected_stimnulated.gexf")
# %%