# %% ---
# IMPORT DATAFRAMES

import pandas as pd
import numpy as np

from scipy.optimize import differential_evolution
from scipy.optimize.optimize import OptimizeResult

import ray

if not ray.is_initialized():
    ray.init()

RANDOM_SEED: int = 23  # Usada para random states
WORKERS: int = 16  # Number of paralel cores by job <= total cores in cluster

# Final agregation proveninete del step04, con multiindices (variation, subsystem, aggregation)
# BUG Degree centrality tiene una cantidad aberrante de NaNs, posiblemente ramas removidas?

agg_df = pd.read_parquet("data/final_aggregation.parquet.gzip")
agg_df.name = "agg_df"
agg_df.sort_index(inplace=True)

# Importa un grupo de reacciones de cuerpos cetonicos. Existen otros grupos tipo branch amminoacids
# la idea de esto es que estas reacciones queden en el ranking mas alto
# osea, son el INPUT Objetivo --> buscamos los pesos para que esta lista sea el top
target_rxns_arr = pd.read_csv("./data/BCAAs_rxns.csv")["."]

# %% ---
# FUNCION A OPTIMIZAR
def get_weighted_ranking(
    weights: list,  # Lista de floats
    agg_df: pd.DataFrame = agg_df,  # Dataframe de valores
    opt_sr: pd.Series = target_rxns_arr,  # Serie a optimizar
    df_loc: tuple = ("NormalizedDelta", "ETC_astrocyte", "quadratic_mean"),
) -> float:
    """Minimizes the ranking for a given series of metabolites"""
    df = agg_df.loc[df_loc, :]
    # Usa el .loc de arriba con un segmento indexado

    # Desempaca los pesos. La cantidad de pesos es igual a las centralidades
    w_vector = w1, w2, w3, w4, w5, w6, w7, w8, w9 = weights

    # Multiplica los pesos por el dataframe, osea cada centralidad por un peso
    C_df_weighted = df * w_vector

    # Promedio por centralidades, osea las centralidades de cada nodo -> hace una combinacion lineal
    C_df_weighted["mean"] = C_df_weighted.mean(axis=1)

    # Rankeado de la optimizacion
    C_df_weighted.sort_values(by=["mean"], inplace=True)  # Ordena la combinacion lineal
    ranking = pd.DataFrame(
        C_df_weighted.index.values
    ).reset_index()  # Sacamos el indice de los nodos
    ranking.columns = ["position", "reaction"]

    # Suma las posiciones de las target reactions
    # por el indice de arriba, estas serian el numero a minimizar
    return ranking[ranking["reaction"].isin(opt_sr)]["position"].values.sum()


# OPTIMIZADOR
@ray.remote
def optimize_weights(
    agg_df: pd.DataFrame = agg_df,
    df_loc: tuple = ("NormalizedDelta", "ETC_astrocyte", "quadratic_mean"),
    b_min: int = 0,
    b_max: int = 1000,
):
    """Runs the optimization for a given set of indexes `df_loc`"""
    np.random.seed(RANDOM_SEED)  # Replicabilidad

    n_centralities: int = agg_df.columns.size  # N of centralities

    # Genera tantos 1 como centralidades existen con
    baseline_rank = get_weighted_ranking([1.0] * n_centralities)

    # Genera una lista de bounds, de [[0,1000],[0,1000]...] 10 veces
    weights_bounds = [(b_min, b_max)] * n_centralities

    # Optimizador inicializado con aleatorios entre [0,1000]
    x0 = np.random.uniform(low=b_min, high=b_max, size=(1, n_centralities))

    # Optimizacion por differential evolution
    # TODO Pasarle el reto de los parametros al optimizador args=()
    solution = differential_evolution(
        get_weighted_ranking,
        args=(agg_df, target_rxns_arr, df_loc),
        seed=RANDOM_SEED,
        x0=x0,
        updating="immediate",
        mutation=(0.07, 1.9),
        bounds=weights_bounds,
        # workers=WORKERS, # BUG Este paralelismo no es compatible con Ray
        workers=1,  # pero con un unico worker funciona?
        maxiter=10 ** 9,
        strategy="rand2exp",
    )

    return [solution, baseline_rank]


# %% ---
# COMPLETE ITERATION

# Get the multiindex minus the 'removed_rxn' level, as unique values
# as a list of tuples used for .loc[index,:]
iterate: list = agg_df.index.droplevel("removed_rxn").unique().tolist()

# TODO #5 Set optimization of centralities as a distribuited computation
optimals = [optimize_weights.remote(agg_df=agg_df, df_loc=t) for t in iterate]
optimals = ray.get(optimals)

# Uses unzip to separate the pair of tuples with results
solutions, baseline_ranks = list(zip(*optimals))

# Dataframe with results; sol.x are the weights for each optimization,
# agg_df.columns are the centralities, and iterate is converted to a
# multiindex, for the first three levels of the aggregated dataframe
daframe = pd.DataFrame(
    [sol.x for sol in solutions],
    columns=agg_df.columns,
    index=pd.MultiIndex.from_tuples(iterate),
)

daframe.index.names = ["Variation", "Subsystem", "Reduction"]

daframe["Ranking improvement"] = list(
    pd.Series(baseline_ranks) - pd.Series([sol.fun for sol in solutions])
)
# %%
# Save the thingi
daframe.to_parquet("data/weights_optimizations.parquet.gzip")

# %%
