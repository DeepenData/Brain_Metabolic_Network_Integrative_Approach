# %% 
from enum import unique
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)

import cobra
import warnings 
import pandas as pd
import networkx as nx
warnings.filterwarnings("ignore")
#Import model from json file
stimulated      = cobra.io.load_json_model(local_working_directory + '/stimulated_2021.json')
# %% 



genes_entrez = [stimulated.genes[g].id for g in range(len(stimulated.genes))]


[r.id for r  in stimulated.genes.get_by_id("1594").reactions]


hola =[]
for g in genes_entrez:
    h = [r.id for r  in stimulated.genes.get_by_id(g).reactions]
    hola.append(h)

  

gene_rxns = pd.DataFrame(zip(hola), genes_entrez, columns=["rxns"])


gene_rxns.to_csv("NAMU_genes_2_rxns_mapping.csv")

# %% 
#Checking bounds in the neuronal glucose uptake rate
stimulated.reactions.get_by_id('GLCt1r').bounds = 3.452112, 3.452112
stimulated.reactions.get_by_id('GLCt1r_Neuron').bounds = 2.120199, 2.120199
stimulated.reactions.get_by_id('DM_atp(c)_Neuron').bounds = 38, 38
stimulated.reactions.get_by_id('ATPS4m').bounds = -1000, -0.01666

#[reaction.id for reaction in stimulated.reactions]
import re
import numpy as np


rxn_list = [reaction.id for reaction in stimulated.reactions]

hh   = [ re.findall(r'(?i)(.*GLC.*|.*DM_atp.*|.*EX_o2.*)', h) for h in rxn_list] 

rxns = [item for sublist in hh for item in sublist]


''
#rxns = ["GLCt1r", "GLCt1r_Neuron", "DM_atp(c)_Neuron"]


bounds = [stimulated.reactions.get_by_id(r).bounds for r in rxns]

my_dict =  dict(zip(rxns , bounds))

pd.DataFrame(my_dict, index=["Lower", "Upper"]).T


rxns_plus = ['NaEX_Neuron','NaKt_Neuron','L-LACt2r_Int','L-LACt2r_Neuron','GLUVESSEC_Neuron',
'GLNtN1_Int','ATPS4m_Neuron','ATPS4m','GLCt1r','GLCt1r_Neuron','ATPtm_Neuron','PYK_Neuron','PYK']

my_rxns = np.unique(rxns+rxns_plus)

#Discrepancies
#print(
#100*abs(1.7 - 1.8823)/1.7, 100*abs(38 - 33.78)/38)



# FBA

#Generate uptake ranges to compute PHPPs
fba_solution                       =  stimulated.optimize()

fba = fba_solution.to_frame()
fba.loc[my_rxns,:]
# %%
