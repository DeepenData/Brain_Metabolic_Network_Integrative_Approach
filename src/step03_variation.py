#!/bin/python
"""Genera un grafo a partir de un modelo de COBRA en formato JSON. Lo exporta en
un formato apto para NetworkX y Graph-tool

INPUTS: 
    "data/centralities.parquet.gzip" : pd.DataFrame
        Centralia
OUTPUT:
    data/NormalizedDelta_df.parquet.gzip : pd.DataFrame
        Centrality variation when comparing unperturbed versus perturbed NAMU
    data/log2Ratio_df.parquet.gzip : pd.DataFrame
        Centrality variation when comparing unperturbed versus perturbed NAMU
"""


#%%
import pandas as pd
import numpy as np
import os

if "src" in os.getcwd():
    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)

# Imports the dataframe with centralities by removed node + baseline
# using a multiindex in level 0 for those.
centralidades_df = pd.read_parquet("data/centralities_long.parquet.gzip")
# check the number of remove reactions

assert len(centralidades_df.index.get_level_values("removed_node").unique()) == 1057  #

# Index sorting
centralidades_df.sort_index(inplace=True)

baseline = centralidades_df.loc["baseline"]

# %% ---
# VARIATIONS - FOLD CHANGE
# is defined as log2(X_initial/X_removed)
def log2RatioCentrality(a_rxn: str) -> pd.DataFrame:

    # SANITY CHECK:  True when all indexes are OK
    assert all(
        baseline.index.values == centralidades_df.loc[a_rxn].index.values
    ), "Index miismatch"

    return np.log2(baseline / centralidades_df.loc[a_rxn])


# holafoldchange = log2RatioCentrality("1a,24,25VITD3Hm")

# VARIATIONS - DELTA
# is defined as (X_initial - X_removed)/X_removed
def NormalizedDeltaCentrality(a_rxn: str) -> pd.DataFrame:
    # SANITY CHECK:  True when all indexes are OK
    assert all(
        baseline.index.values == centralidades_df.loc[a_rxn].index.values
    ), "Index miismatch"

    return (baseline - centralidades_df.loc[a_rxn]) / baseline.loc[a_rxn]


def get_centrality_variation(variation_function, rxns: list) -> pd.DataFrame:

    if "baseline" in rxns:
        rxns.remove("baseline")

    my_iter = map(variation_function, rxns)
    Centralities_list = []
    for i_df, rxn in zip(my_iter, all_rxns):

        i_df["removed_rxn"] = rxn

        i_df.reset_index(inplace=True)
        i_df.set_index(["removed_rxn", "metabolite"], inplace=True)
        Centralities_list.append(i_df)

    return pd.concat(Centralities_list, axis=0)


all_rxns = list(
    np.unique(centralidades_df.index.get_level_values("removed_node").values)
)
#%%
# Centrality variation when comparing unperturbed versus perturbed NAMU
NormalizedDelta_df = get_centrality_variation(NormalizedDeltaCentrality, all_rxns)
NormalizedDelta_df.to_parquet("data/NormalizedDelta_df.parquet.gzip")

#
log2Ratio_df = get_centrality_variation(log2RatioCentrality, all_rxns)
log2Ratio_df.to_parquet("data/log2Ratio_df.parquet.gzip")

# %%
