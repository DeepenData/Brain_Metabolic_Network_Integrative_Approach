# %% 
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
#bounds in the neuronal glucose uptake rate
stimulated.reactions.get_by_id('GLCt1r').bounds =3.45501146, 3.45501146 
stimulated.reactions.get_by_id('GLCt1r_Neuron').bounds = 2.1177, 2.1177
stimulated.reactions.get_by_id('DM_atp(c)_Neuron').bounds = 38, 38
stimulated.reactions.get_by_id('ATPS4m').bounds = -1000, -0.01666

#Generate uptake ranges to compute PHPPs
fba_solution                       =  stimulated.optimize()
my_rxns = ['GLCt1r_Neuron', 'GLCt1r','ATPS4m','DM_atp(c)_Neuron','NaEX_Neuron',
'NaKt_Neuron','L-LACt2r_Int','L-LACt2r_Neuron',
'GLUVESSEC_Neuron','GLNtN1_Int','ATPS4m_Neuron']
fba = pd.DataFrame(
[(stimulated.reactions.get_by_id(r).id, 
stimulated.reactions.get_by_id(r).name,
stimulated.reactions.get_by_id(r).reaction,
stimulated.reactions.get_by_id(r).bounds,
fba_solution.to_frame().loc[r,].fluxes)  for r in my_rxns])
fba.index = list(range(1,12))
fba.columns = ["ID","Name","Reaction","Bounds","Flux"]
fba
# %% 
optimal_astrocytic_glucose_uptake, optimal_overall_oxygen_uptake_rate  =  fba_solution["GLCt1r"], -1*fba_solution["EX_o2(e)"]

def generate_uptake_range(optimal, first_step, last_step, direction, intervals=100):
    import numpy as np
    lower   = np.linspace(first_step, optimal, num=intervals, endpoint=True)
    upper   = np.linspace(optimal, last_step*optimal,num=intervals, endpoint=True)
    myrange = direction*np.unique(np.concatenate(( lower, upper)))
    return myrange 


range_of_glucose_uptake_rates = generate_uptake_range(optimal_astrocytic_glucose_uptake, 1, 1.7,1) 
range_of_oxygen_uptake_rates  = generate_uptake_range(optimal_overall_oxygen_uptake_rate, 13, 1.5, -1)
# %%
import numpy as np
#Create empty objects
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

#Compute PHPP
for i in range(0,l):
    for j in range(0,l):
        stimulated.reactions.get_by_id("EX_o2(e)").bounds = (range_of_oxygen_uptake_rates[i], 0)
        stimulated.reactions.get_by_id("GLCt1r").bounds   = (range_of_glucose_uptake_rates[j],range_of_glucose_uptake_rates[j])        
        sol              = stimulated.optimize()  
        objective[i,j]        = sol.objective_value
        NaEX_Neuron[i,j]      =sol.fluxes['NaEX_Neuron'] 
        NaKt_Neuron[i,j]      =sol.fluxes['NaKt_Neuron']
        L_LACt2r_Int[i,j]     =sol.fluxes['L-LACt2r_Int']
        L_LACt2r_Neuron[i,j]  =sol.fluxes['L-LACt2r_Neuron']
        GLUVESSEC_Neuron[i,j] =sol.fluxes['GLUVESSEC_Neuron']
        GLNtN1_Int[i,j]       =sol.fluxes['GLNtN1_Int']
        ATPS4m_Neuron[i,j]    =sol.fluxes['ATPS4m_Neuron']
        ATPS4m[i,j]           =sol.fluxes['ATPS4m']
        glucose_uptake[i,j]           =sol.fluxes['GLCt1r']
        GLCt1r_Neuron[i,j]    =sol.fluxes['GLCt1r_Neuron']
        ATPtm_Neuron[i,j]     =sol.fluxes['ATPtm_Neuron']
        PYK_Neuron[i,j]       =sol.fluxes['PYK_Neuron']
        PYK[i,j]              =sol.fluxes['PYK']
        #EX_o2_e[i,j]          =sol.fluxes[rxns.index('EX_o2(e)')]        
        #EX_glc_e[i,j]          =sol.fluxes[rxns.index('EX_glc(e)')]        
        oxygen_uptake[i,j] = stimulated.reactions.get_by_id("EX_o2(e)").bounds[0]

# %% 
figure_path = local_working_directory + "/Figure_3/"
#Save outputs
from numpy import savetxt
savetxt(figure_path+'glucose_uptake.csv', glucose_uptake, delimiter=',')
savetxt(figure_path+'oxygen_uptake.csv', -1*oxygen_uptake, delimiter=',')
#savetxt(figure_path+'GLNtN1_Int.csv', GLNtN1_Int, delimiter=',')
savetxt(figure_path+'GLUVESSEC_Neuron.csv', GLUVESSEC_Neuron, delimiter=',')
savetxt(figure_path+'L_LACt2r_Neuron.csv', L_LACt2r_Neuron, delimiter=',')
savetxt(figure_path+'L_LACt2r_Int.csv', -1*L_LACt2r_Int, delimiter=',')
#savetxt(figure_path+'NaKt_Neuron.csv', NaKt_Neuron, delimiter=',')
savetxt(figure_path+'NaEX_Neuron.csv', NaEX_Neuron, delimiter=',')
savetxt(figure_path+'objective.csv', objective, delimiter=',')
#Compute yield
Y_o2_glc    = -1*oxygen_uptake/(GLCt1r_Neuron+glucose_uptake)

Y_ATP_glc   =  (ATPtm_Neuron+PYK_Neuron)/(GLCt1r_Neuron+glucose_uptake)

savetxt(figure_path+'Y_o2_glc.csv', Y_o2_glc, delimiter=',')
savetxt(figure_path+'Y_ATP_glc.csv', Y_ATP_glc, delimiter=',')
# %% Make a dataframe with the FBA results
fba_solution_df = fba_solution.to_frame()
rxns = [stimulated.reactions[i].reaction for i in range(0,len(stimulated.reactions))]
names = [stimulated.reactions[i].name for i in range(0,len(stimulated.reactions))]
FBA_results = pd.DataFrame(list(zip(fba_solution_df.index, names ,rxns, fba_solution_df.fluxes, fba_solution_df.reduced_costs)), 
columns= ['ID', 'Name','Reaction','Flux','Sensitivity'])
FBA_results.to_csv(figure_path+"FBA_results.csv")