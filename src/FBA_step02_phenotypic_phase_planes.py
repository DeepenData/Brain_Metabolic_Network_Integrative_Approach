#%%
#import pandas as pd
import modin.pandas as pd
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
#SETTING BOUNDS from experimental fluxes
stimulated.reactions.get_by_id('GLCt1r').bounds = 3.452112, 3.452112
stimulated.reactions.get_by_id('GLCt1r_Neuron').bounds = 2.120199, 2.120199
stimulated.reactions.get_by_id('DM_atp(c)_Neuron').bounds = 38, 38
stimulated.reactions.get_by_id('ATPS4m').bounds = -1000, -0.01666
#CHECK FBA solution for neurotransmission-related reactions
rxn_list      = [reaction.id for reaction in stimulated.reactions]
rxn_name_list = [reaction.name for reaction in stimulated.reactions]
rxn_stoich_list = [reaction.reaction for reaction in stimulated.reactions]
# %%
FBA_step01_neurotransmission_related_fluxes = pd.read_parquet("data/FBA_step01_neurotransmission_related_fluxes.parquet.gzip")


optimal_astrocytic_glucose_uptake, optimal_overall_oxygen_uptake_rate  =  FBA_step01_neurotransmission_related_fluxes.loc['GLCt1r']['fluxes'], -1*FBA_step01_neurotransmission_related_fluxes.loc['EX_o2(e)']['fluxes']

def get_uptake_range(optimal, first_step, last_step, direction, intervals=100):
    import numpy as np
    lower   = np.linspace(first_step, optimal, num=intervals, endpoint=True)
    upper   = np.linspace(optimal, last_step*optimal,num=intervals, endpoint=True)
    return direction*np.unique(np.concatenate(( lower, upper))) 

range_of_glucose_uptake_rates = get_uptake_range(optimal_astrocytic_glucose_uptake, 1, 1.7,    1)
range_of_oxygen_uptake_rates  = get_uptake_range(optimal_overall_oxygen_uptake_rate, 13, 1.5, -1)

l                   = len(range_of_glucose_uptake_rates)
NaEX_Neuron         =  np.zeros((l,l))
NaKt_Neuron         =  np.zeros((l,l))
L_LACt2r_Int        =  np.zeros((l,l))
L_LACt2r_Neuron     =  np.zeros((l,l))
GLUVESSEC_Neuron    =  np.zeros((l,l))
GLNtN1_Int          =  np.zeros((l,l))
ATPS4m_Neuron       =  np.zeros((l,l))
ATPS4m              =  np.zeros((l,l))
glucose_uptake      =  np.zeros((l,l))
GLCt1r_Neuron      =  np.zeros((l,l))
#EX_o2_e            =  np.zeros((l,l)) 
#EX_glc_e           =  np.zeros((l,l)) 
oxygen_uptake       =  np.zeros((l,l))
objective           =  np.zeros((l,l)) 

ATPtm_Neuron =  np.zeros((l,l))
PYK_Neuron   =  np.zeros((l,l))
PYK          =  np.zeros((l,l)) 

def compute_fba(i,j, stimulated):
    stimulated.reactions.get_by_id("EX_o2(e)").bounds = (range_of_oxygen_uptake_rates[i], 0)
    stimulated.reactions.get_by_id("GLCt1r").bounds   = (range_of_glucose_uptake_rates[j],range_of_glucose_uptake_rates[j])        
    sol                   = stimulated.optimize()  
    objective[i,j]        = sol.objective_value
    NaEX_Neuron[i,j]      =sol.fluxes['NaEX_Neuron'] 
    NaKt_Neuron[i,j]      =sol.fluxes['NaKt_Neuron']
    L_LACt2r_Int[i,j]     =sol.fluxes['L-LACt2r_Int']
    L_LACt2r_Neuron[i,j]  =sol.fluxes['L-LACt2r_Neuron']
    GLUVESSEC_Neuron[i,j] =sol.fluxes['GLUVESSEC_Neuron']
    GLNtN1_Int[i,j]       =sol.fluxes['GLNtN1_Int']
    ATPS4m_Neuron[i,j]    =sol.fluxes['ATPS4m_Neuron']
    ATPS4m[i,j]           =sol.fluxes['ATPS4m']
    glucose_uptake[i,j]    =sol.fluxes['GLCt1r']
    GLCt1r_Neuron[i,j]    =sol.fluxes['GLCt1r_Neuron']
    ATPtm_Neuron[i,j]     =sol.fluxes['ATPtm_Neuron']
    PYK_Neuron[i,j]       =sol.fluxes['PYK_Neuron']
    PYK[i,j]              =sol.fluxes['PYK']
    #EX_o2_e[i,j]          =sol.fluxes[rxns.index('EX_o2(e)')]
    #EX_glc_e[i,j]          =sol.fluxes[rxns.index('EX_glc(e)')]        
    oxygen_uptake[i,j] = stimulated.reactions.get_by_id("EX_o2(e)").bounds[0]
    return objective, NaEX_Neuron, NaKt_Neuron, L_LACt2r_Int, L_LACt2r_Neuron, GLUVESSEC_Neuron, GLNtN1_Int, ATPS4m_Neuron, ATPS4m, glucose_uptake,  \
           GLCt1r_Neuron,  ATPtm_Neuron, PYK_Neuron,  PYK, oxygen_uptake
 
# %%
#Compute PHPP
for i in range(l):
    for j in range(l):
        objective, NaEX_Neuron, NaKt_Neuron, L_LACt2r_Int, L_LACt2r_Neuron, GLUVESSEC_Neuron, GLNtN1_Int, \
        ATPS4m_Neuron, ATPS4m, glucose_uptake, GLCt1r_Neuron,  ATPtm_Neuron, PYK_Neuron,  PYK, oxygen_uptake = compute_fba(i,j, stimulated)
        
def get_df(arr, name):
    
    df          = pd.DataFrame(arr)     
    df['index'] = name
    df.name     = name
    return df


phpps       = [objective, NaEX_Neuron, NaKt_Neuron, L_LACt2r_Int, L_LACt2r_Neuron, GLUVESSEC_Neuron, GLNtN1_Int, \
                ATPS4m_Neuron, ATPS4m, glucose_uptake, GLCt1r_Neuron,  ATPtm_Neuron, PYK_Neuron,  PYK, oxygen_uptake] 

phpps_index = ['objective', 'NaEX_Neuron', 'NaKt_Neuron', 'L_LACt2r_Int', 'L_LACt2r_Neuron', 'GLUVESSEC_Neuron', 'GLNtN1_Int', \
                 'ATPS4m_Neuron', 'ATPS4m', 'glucose_uptake', 'GLCt1r_Neuron',  'ATPtm_Neuron', 'PYK_Neuron',  'PYK', 'oxygen_uptake'] 

phpp_iter         = map(get_df, phpps, phpps_index)
all_phpps         =  pd.concat(list(phpp_iter), axis=0).set_index('index')
all_phpps.columns =  [str(i) for i in all_phpps.columns.values]
all_phpps.to_parquet("data/FBA_step02_all_phpps.parquet.gzip")
# %%