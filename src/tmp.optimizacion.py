
#%%
import pandas as pd
from scipy.optimize import differential_evolution
import os
from itertools import product
import numpy as np
my_cpus = os.cpu_count()
centralities   = ['bet', 'clo', 'page','bet2', 'clo2', 'page2']
n_reactions    = 11000
n_centralities = len(centralities)
arr1           = 'r'
arr2           = range(n_reactions)
names_iter     = product(arr1, arr2)
rxns           = [''.join(map(str, v)) for v in names_iter]

def make_random_data(cols=n_centralities, rows=n_reactions):
    while True:
        yield np.random.uniform(low=.01, high=2, size=(rows, cols))       #crea una distribución random uniforme

rand_gen = make_random_data()

C_df = pd.DataFrame(next(rand_gen), index = rxns, columns=centralities)

def objFun(w):
    w_vector = w1, w2, w3, w4, w5, w6 = w

    C_df_weighted         = C_df*np.array(w_vector)
    C_df_weighted['mean'] = C_df_weighted.mean(axis=1)
    C_df_weighted.sort_values(by=['mean'], inplace=True)
    ranking         = pd.DataFrame(C_df_weighted.index.values).reset_index()
    ranking.columns = ['position', 'reaction']
    score           = ranking[ranking.reaction.isin(target_rxns)].position.values.sum()

    return score

r_min, r_max = 0, 100
bounds = [[r_min, r_max]]
bounds_x = [] # 
n_vars = n_centralities
for _ in range(n_vars):
    bounds_x.extend(bounds)

# %%
x0 =  np.random.uniform(low=r_min, high=r_max, size=(1, n_centralities))  
target_rxns = ['r7', 'r3', 'r0', 'r88', 'r100','r23']

solution = differential_evolution(
        objFun,
        x0= x0, updating='immediate',mutation=(0.1, 1.8), 
        bounds = bounds_x, workers=4*my_cpus,
        maxiter=10**9, strategy = 'rand2exp')      

solution


# %% testing
import numpy as np

np.array([1,2,3,4]).astype(int)




'''w           = [54.02691898, 36.99043141, 30.03487773]
target_rxns = ['r7', 'r3', 'r0']

C_df_weighted         = C_df*np.array(w)
C_df_weighted['mean'] = C_df_weighted.mean(axis=1)
C_df_weighted.sort_values(by=['mean'], inplace=True)

ranking         = pd.DataFrame(C_df_weighted.index.values).reset_index()
ranking.columns = ['position', 'reaction']

A =  ranking[ranking.reaction.isin(target_rxns)]
score           = ranking[ranking.reaction.isin(target_rxns)].position.values.sum()
print('score',score, '\n',A)'''