#!/bin/python
"""Genera un grafo a partir de un modelo de COBRA en formato JSON. Lo exporta en
un formato apto para NetworkX y Graph-tool

INPUTS: 
    graphml_path : (str) Direccion del .graphml, eg. "./tmp/graph.graphml"
    nodos_remover : (list) Lista de str de nodos a remover del grafo independientemente. Puede ser lista de listas
OUTPUT: 
    data/centralities.parquet.gzip : (DataFrame) DataFrame de Pandas con centralidades, y multiindice baseline + removidos
"""
#%%
import pandas as pd
from pandas.core.frame import DataFrame
import numpy as np
import os

if "src" in os.getcwd():
    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)
    
NormalizedDelta_df = pd.read_parquet("data/NormalizedDelta_df.parquet.gzip")
NormalizedDelta_df.sort_index(inplace=True)
NormalizedDelta_df.name = "NormalizedDelta"
NormalizedDelta_df.replace([np.inf, -np.inf], np.nan, inplace=True)

log2Ratio_df = pd.read_parquet("data/log2Ratio_df.parquet.gzip")
log2Ratio_df.sort_index(inplace=True)
log2Ratio_df.name = "log2Ratio"
log2Ratio_df.replace([np.inf, -np.inf], np.nan, inplace=True)

FBA_sensitive = pd.read_csv('data/FBA_sensitivities.csv')


objective = [
"GLNtN1_Int", 
"GLUVESSEC_Neuron", 
"L-LACt2r_Int",
"L-LACt2r_Neuron",
"NaKt_Neuron"]


# TODO Subsystems should be given from a config file, not hardcoded
subsystems = {"Objective": objective, 
    'FBA_sensitive': FBA_sensitive.iloc[:,0].values,
    "Glycolysis_astrocyte": [
        "PGM",
        "ACYP",
        "PGI",
        "PGK",
        "PYK",
        "HEX1",
        "DPGase",
        "TPI",
        "PFK",
        "ENO",
        "GAPD",
        "DPGM",
        "FBA",
        "G3PD2m",
    ],
    "Glycolysis_neuron": [
        "ACYP_Neuron",
        "DPGM_Neuron",
        "DPGase_Neuron",
        "ENO_Neuron",
        "FBA_Neuron",
        "G3PD2m_Neuron",
        "GAPD_Neuron",
        "HEX1_Neuron",
        "PFK_Neuron",
        "PGI_Neuron",
        "PGK_Neuron",
        "PGM_Neuron",
        "PYK_Neuron",
        "TPI_Neuron",
    ],
    "ETC_neuron": [
        "ATPS4m_Neuron",
        "CYOOm2_Neuron",
        "CYOR-u10m_Neuron",
        "NADH2-u10m_Neuron",
        "PPA_Neuron",
        "PPAm_Neuron",
    ],
    "ETC_astrocyte": [
        "PPAm",
        "ATPS4m",
        "CYOOm2",
        "CYOR-u10m",
        "NADH2-u10m",
        "PPA",
    ],
}
import itertools

key_to_extract = {"Glycolysis_neuron","Glycolysis_astrocyte", "ETC_astrocyte", "ETC_neuron"}
original_list  = [value for key,value in subsystems.items() if key in key_to_extract]


energy_metabolism = list(itertools.chain(*original_list))


subsystems.update({"energy_metabolism":energy_metabolism })

# log2Ratio_df[log2Ratio_df.index.isin(subsystems["Glycolysis_astrocyte"], level=1)]

# Programar una funcion que tome un Df de forma level_0, level_1...
# [removed_rxn, metabolite] x [harmonic_centrality  eigenvector_centrality  closeness_centrality  information_centrality]
# y haga un promedio por level 1, devolviendo un dataframe
# [removed_rxn] x [harmonic_centrality  eigenvector_centrality  closeness_centrality  information_centrality]



# %% --- DEFINE COSAS DE PROMEDIOS

arr = np.array([0,4,4,np.nan,0,1])
arr[np.nonzero(arr)]

# %% 
import numpy as np
from  numba import jit
from  scipy.stats import gmean, hmean
from  numpy import nanmean

# Define una función helper para cuadraticos (Root Square Mean). No en scipy base
# @jit(nopython=True)
#BUG division por cero en harmonic_mean y geometric_mean


def quadratic_mean(array):
    # array = array.values
    return np.sqrt(np.nanmean(array * array))

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
# @jit(nopython=True)
def geometric_mean(array):
    array = array[~np.isnan(array)] # Filtra los NaN
    
    if  len(array) == 0:
        return np.nan
    
    fcorr = len(array[array > 0])/len(array) # Factor de corrección
    array = array[array > 0] # Selecciona solo mayores a cero 

    
    try: 
        gmean = ( np.prod( array )**(1 / len(array)) ) * fcorr # Geometrica por factor de corrección
    except: 
        return np.nan
    else : 
        return gmean

# Requiere definir un helper con un factor de correction NoCeros/ValoresTotales
# @jit(nopython=True)
def harmonic_mean(array):
    array = array[~np.isnan(array)]  # Filtra los NaN
    if  len(array) == 0:
        return np.nan
    fcorr = len(array[array > 0]) / len(array)  # Factor de corrección
    array = array[array > 0]  # Selecciona solo mayores a cero
    
    if  all(array == 0):
        return 0
    try:
        hmean = (
            len(array) / np.sum(1 / array) * fcorr
        )  # Geometrica por factor de corrección
    except:
        return np.nan
    else:
        return hmean

# TODO Mover esto a un test_step04_...
def gmean_scipy(a_series):
    return gmean(a_series[~a_series.isna()])

def hmean_scipy(a_series):
    return hmean(a_series[~a_series.isna()])


def nanmean_numpy(a_series):
    return nanmean(a_series[~a_series.isna()])


# Esto no es parte del test!
#geometric = NormalizedDelta_df.groupby("removed_rxn").aggregate(gmean_scipy)

# log2Ratio_df[log2Ratio_df.index.isin(subsystems["Glycolysis_astrocyte"], level=1)]


def make_aggregation_by_subsystem(dataframe: pd.DataFrame, aggregation, subsystem: str):
    """
    [summary]

    [extended_summary]

    Parameters
    ----------
    dataframe : pd.DataFrame
        [description]
    aggregation : function
        Operator, for computation over the "removed_rxn"

    Returns
    -------
    pd.DataFrame
        [description]
    """
    df_name = dataframe.name
    dataframe = dataframe[dataframe.index.isin(subsystems[subsystem], level=1)]
    dataframe = dataframe.groupby("removed_rxn").aggregate(aggregation)
    dataframe.reset_index(inplace=True)
    dataframe["aggregation"] = aggregation.__name__
    dataframe["subsystem"] = subsystem
    dataframe["variation"] = df_name
    dataframe.set_index(
        ["variation", "subsystem", "aggregation", "removed_rxn"], inplace=True
    )
    dataframe.sort_index(inplace=True)
    return dataframe


make_aggregation_by_subsystem(
    NormalizedDelta_df, quadratic_mean, "Glycolysis_astrocyte"
)


import itertools as itt


all_subsystems = list(subsystems.keys())
#%%
#
args_iterator = itt.product(
    [log2Ratio_df, NormalizedDelta_df],
    [quadratic_mean, gmean_scipy, nanmean_numpy],#gmean_scipy, nanmean_numpy, quadratic_mean 
    all_subsystems,
)

agg_list = []
for df, func, subs in args_iterator:
    tmp = make_aggregation_by_subsystem(df, func, subs)
    agg_list.append(tmp)

final_aggregation = pd.concat(agg_list, axis=0)
#%%
final_aggregation.to_parquet("data/final_aggregation.parquet.gzip")

# %%
