# NAMU2021b
## NAMU2023a

1. COBRA model -> Networkx ?
2. Networkx -> Centralities ?

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
