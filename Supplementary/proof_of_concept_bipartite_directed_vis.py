# %% --- 
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
supp_path    = local_working_directory + "/Supplementary/"
supp_path_results    = local_working_directory + "/Supplementary/" + "results/"
from aa_functions import *
import numpy as np
from cobra import *
from   cobra.util.array import create_stoichiometric_matrix
import pandas as pd
# %% --- 
model = create_toy_metabolism()
S_matrix = create_stoichiometric_matrix(model)

pd.DataFrame(S_matrix)
# %%
