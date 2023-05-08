
# Source documentation

This is mainly written in Python, with some R code for generation beautifull ggplot figures. 

The `conda_NAMU2021b.yaml` file has the current environment and installed packages. 

## Workflow

<!--This should be implemented in Snakemake-->

- `step_01`: converts the model from a `.json` to a `.graphml`, representing a proyected graph of the reactions, with edges corresponding to metabolites. 

- `step_02`: computes centralities using NetworkX and Graph-tool, the later being prefered as for the faster C implementation. 
  Nodes are removed as to simulate a reaction knock-out, akin to a loss-of-function mutation. 

- `step_03`: computes the fold-change and deltas in centrality for the removed nodes. 

- `step_04`: makes the ggplot figures. 
