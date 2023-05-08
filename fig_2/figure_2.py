#%%
import pandas as pd
import sklearn as sk
from sklearn.cluster import *


if "fig_1" in os.getcwd():
    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)
#%%
def make_clustering(df, quantile = .2):
    df_copy = df.copy()
    X                             = df_copy.AbsoluteOP.values.reshape(-1, 1)
    bandwidth                     = estimate_bandwidth(X, quantile=quantile, n_samples=1000, random_state=0, n_jobs=-1)
    mean_shift_clusters           = mean_shift(X, bandwidth = bandwidth, seeds=None, bin_seeding=False, 
                                    min_bin_freq=1, cluster_all=True, max_iter=500, n_jobs=-1)
    _, df_copy['mean_shift_clusters']  = mean_shift_clusters
    return(df_copy)

def get_intervals(cluster, clusterized):
    
    A = clusterized.loc[clusterized.mean_shift_clusters == cluster].AbsoluteOP
    return (cluster,A.min(),A.shape[0],A.max())

def evaluate_clustering(data, quantile):
    clusterized = make_clustering(data, quantile)
    clusters    =  clusterized.mean_shift_clusters.unique()

    return clusterized, [get_intervals(c, clusterized) for c in clusters]
  #%%
clusterized_Astrocyte, intervals_Astrocyte = \
evaluate_clustering(pd.read_csv("data/Astrocyte_Absolute_Optimality.csv"), quantile=.375 )

clusterized_Neuron, intervals_Neuron      = \
evaluate_clustering(pd.read_csv("data/Neuron_Absolute_Optimality.csv"), quantile= .375)

clusterized_Astrocyte.to_csv("data/Astrocyte_Absolute_Optimality_clusterized.csv")
clusterized_Neuron.to_csv("data/Neuron_Absolute_Optimality_clusterized.csv")

intervals_Astrocyte
intervals_Neuron
#%%
import pandas as pd
import sklearn as sk
from sklearn.cluster import *
import os

if "fig_1" in os.getcwd():
    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)

def make_clustering(df, quantile = .2):
    df_copy = df.copy()
    X                             = df_copy.AbsoluteC.values.reshape(-1, 1)
    bandwidth                     = estimate_bandwidth(X, quantile=quantile, n_samples=1000, random_state=123, n_jobs=-1)
    mean_shift_clusters           = mean_shift(X, bandwidth = bandwidth, seeds=None, bin_seeding=False, 
                                    min_bin_freq=1, cluster_all=True, max_iter=500, n_jobs=-1)
    _, df_copy['mean_shift_clusters']  = mean_shift_clusters
    return(df_copy)

def get_intervals(cluster, clusterized):
    
    A = clusterized.loc[clusterized.mean_shift_clusters == cluster].AbsoluteC
    return (cluster,A.min(),A.shape[0],A.max())

def evaluate_clustering(data, quantile):
    clusterized = make_clustering(data, quantile)
    clusters    =  clusterized.mean_shift_clusters.unique()

    return clusterized, [get_intervals(c, clusterized) for c in clusters]


clusterized_Astrocyte, intervals_Astrocyte = \
evaluate_clustering(pd.read_csv("data/Astrocyte_Absolute_Centrality.csv"), quantile=.18 )

clusterized_Neuron, intervals_Neuron      = \
evaluate_clustering(pd.read_csv("data/Neuron_Absolute_Centrality.csv"), quantile= .18)

clusterized_Astrocyte.to_csv("data/Astrocyte_Absolute_Centrality_clusterized.csv")
clusterized_Neuron.to_csv("data/Neuron_Absolute_Centrality_clusterized.csv")



intervals_Astrocyte
intervals_Neuron
#%%
clusterized_Astrocyte, intervals_Astrocyte = \
evaluate_clustering(pd.read_csv("data/Astrocyte_Absolute_Centrality_Sensitive.csv"), quantile=.274 )

clusterized_Neuron, intervals_Neuron      = \
evaluate_clustering(pd.read_csv("data/Neuron_Absolute_Centrality_Sensitive.csv"), quantile= .274)



intervals_Astrocyte

intervals_Neuron


#%%
clusterized_Astrocyte.to_csv("data/Astrocyte_Absolute_Centrality_Sensitive_clusterized.csv")
clusterized_Neuron.to_csv("data/Neuron_Absolute_Centrality_Sensitive_clusterized.csv")
