# Metabolic switch in the aging astrocyte supported via integrative approach comprising network and transcriptome analyses 
## Abstract

Dysregulated central-energy metabolism is a hallmark of brain aging. Supplying enough energy for neurotransmission relies on the neuron-astrocyte metabolic network. To identify genes contributing to age-associated brain functional decline, we formulated an approach to analyze the metabolic network by integrating flux, network structure and transcriptomic databases of neurotransmission and aging. Our findings support that during brain aging: (1) The astrocyte undergoes a metabolic switch from aerobic glycolysis to oxidative phosphorylation, decreasing lactate supply to the neuron, while the neuron suffers intrinsic energetic deficit by downregulation of Krebs cycle genes, including mdh1 and mdh2 (Malate-Aspartate Shuttle); (2) Branched-chain amino acid degradation genes were downregulated, identifying dld as a central regulator; (3) Ketone body synthesis increases in the neuron, while the astrocyte increases their utilization, in line with neuronal energy deficit in favor of astrocytes. We identified candidates for preclinical studies targeting energy metabolism to prevent age-associated cognitive decline.

Keywords: astrocyte; brain aging; flux balance analysis; network centrality; neuron.


## Direct Deps

```
- pandas
- numpy
- sklearn
- cobra
- networkx
- matplotlib
- graph_tool
- ray
- numba
```

Desde `podman pull docker.io/library/python`, 
usamos la [imagen oficial de Python](https://hub.docker.com/_/python/). 

## Workflow

1. Read the `JSON` model (from `data/models`) to a `.graphml`
2. The proyected graph is saved in `data/graphs`
3. Baseline centrality, and perturbated centralities computation, using `src/centrality.py`
4. Compute the fold-changes, deltas, for several aritmetic, geometric means. And, the same but for ranked.
5. Rankings for _keto bodies_ and _BCA_
    1. Optimize weight of each centrality, thus maximizing _Keto bodies_, _BCA_, and _both_ subsystems. <!--TODO: BAD EXPLANATION-->
<!--TODO: Redo the analysis using a weighted (via FBA) network. 
6. From 1, redo using a weighted network using the FBA analysis 
    1. Arches is the same number of metabolites  
    2. Use named arches, as the same metabolites  
-->

<!--TODO: compute this centralities
1. Eigenvector centrality
2. Katz centrality
    Alt, Pagerank
3. Closeness centrality
    Alt, variaciones que acepten peso
4. Information centrality
5. Betweenes centrality
6. Comunicability (only NX)
-->
