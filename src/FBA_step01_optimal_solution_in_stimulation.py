#%%
import pandas as pd
import numpy as np
import os
import cobra
import warnings 
import re
import numpy as np
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
warnings.filterwarnings("ignore")
if "src" in os.getcwd():
    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)    
#load model
stimulated      = cobra.io.load_json_model('data/models/stimulated_NAMU.json')    
#%%


#SETTING BOUNDS from experimental fluxes
stimulated.reactions.get_by_id('GLCt1r').bounds = 3.452112, 3.452112
stimulated.reactions.get_by_id('GLCt1r_Neuron').bounds = 2.120199, 2.120199
stimulated.reactions.get_by_id('DM_atp(c)_Neuron').bounds = 38, 38
stimulated.reactions.get_by_id('ATPS4m').bounds = -1000, -0.01666
#CHECK FBA solution for neurotransmission-related reactions
rxn_list      = [reaction.id for reaction in stimulated.reactions]
objective = [
"GLNtN1_Int", 
"GLUVESSEC_Neuron", 
"L-LACt2r_Int",
"L-LACt2r_Neuron",
"NaKt_Neuron"]
set(objective).issubset(set(rxn_list))


#%%
rxn_name_list = [reaction.name for reaction in stimulated.reactions]
rxn_stoich_list = [reaction.reaction for reaction in stimulated.reactions]
hh       = [re.findall(r'(?i)(.*GLCt1r.*|.*DM_atp.c.*Neuron.*)', h) for h in rxn_list] 
rxns     = [item for sublist in hh for item in sublist]
bounds   = [stimulated.reactions.get_by_id(r).bounds for r in rxns]
my_dict  =  dict(zip(rxns , bounds))
#pd.DataFrame(my_dict, index=["Lower", "Upper"]).T
rxns_plus                        = ['EX_o2(e)','NaEX_Neuron','NaKt_Neuron','L-LACt2r_Int','L-LACt2r_Neuron','GLUVESSEC_Neuron','GLNtN1_Int','ATPS4m_Neuron','ATPS4m','GLCt1r','GLCt1r_Neuron','ATPtm_Neuron','PYK_Neuron','PYK']
my_rxns                          = np.unique(rxns+rxns_plus)
fba_solution_df                  =  stimulated.optimize().to_frame()
fba_solution_df['name']          = rxn_name_list
fba_solution_df['stoichiometry'] = rxn_stoich_list

fba_solution_df.filter(regex='^(?!DM_|sink|EX_)', axis=0).to_parquet("data/FBA_step01_FBA_solution.parquet.gzip")
neurotransmission_related_fluxes = fba_solution_df.loc[my_rxns,:]
neurotransmission_related_fluxes.to_parquet("data/FBA_step01_neurotransmission_related_fluxes.parquet.gzip")
sensitivities = fba_solution_df.filter(regex='^(?!DM_|sink|EX_)', axis=0).sort_values('reduced_costs', 
                                                                                      ascending=False).query('abs(reduced_costs)>1e-14')


# %%
sensitivities.to_parquet("data/FBA_step01_non_zero_sensitivities.parquet.gzip")
# %%
