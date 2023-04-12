# %% 
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
supp_path_results    = local_working_directory + "/Supplementary/" + "results/"

import cobra
import warnings 
import pandas as pd
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
import dataframe_image as dfi
dfi.export(fba, supp_path_results + 'namu_fba.png')