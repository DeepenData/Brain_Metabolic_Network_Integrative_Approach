setwd("~/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism/Supplementary")
source('R_functions.R')
library(tidyverse)
effects_on_sub0 <-  get_results("./results/proj_two_core_diRandom_0.xlsx", "fold_change_geometric","./results/proj_two_core_diRandom_0.txt", "./results/proj_two_core_diRandom_1.txt" )
effects_on_sub1 <-  get_results("./results/proj_two_core_diRandom_1.xlsx", "fold_change_geometric","./results/proj_two_core_diRandom_0.txt", "./results/proj_two_core_diRandom_1.txt" )

effects_on_sub0[[1]] %>% colnames() -> Ns
colnames(effects_on_sub0[[1]]) <- paste0(Ns, '_on_sub0')
effects_on_sub1[[1]] %>% colnames() -> Ns
colnames(effects_on_sub1[[1]]) <- paste0(Ns, '_on_sub1')

get_nodes_by_subgraph(effects_on_sub0[[1]], effects_on_sub0[[2]]) -> sub_0_nodes_by_subgraph
get_nodes_by_subgraph(effects_on_sub1[[1]], effects_on_sub1[[2]]) -> sub_1_nodes_by_subgraph

rbind(sub_0_nodes_by_subgraph$subgraph_0, sub_0_nodes_by_subgraph$subgraph_1) %>% as.data.frame() %>% 
  rename(subgraph_from_sub0 = subgraph) %>% rownames_to_column("node") -> effects_on_sub0_with_types

rbind(sub_1_nodes_by_subgraph$subgraph_0, sub_1_nodes_by_subgraph$subgraph_1) %>% as.data.frame()%>% 
  rename(subgraph_from_sub1 = subgraph)  %>% rownames_to_column("node")-> effects_on_sub1_with_types

first_join <- inner_join(effects_on_sub0_with_types, effects_on_sub1_with_types) 
first_join %>% column_to_rownames("node") %>% select_if(is.character) %>% .[[1]] -> node_types
first_join %>% column_to_rownames("node") %>% select_if(is.double)-> all_fc_centralities
rowSums(abs(all_fc_centralities)) -> total.abs.centrality



library(WGCNA)
all_fc_centralities %>% transposeBigData   %>% cor(method = "kendall") -> corrM
library(ComplexHeatmap)
library(magrittr)
node_types %<>% str_replace_all(".*0.*", "Nodes from yellow subgraph")
node_types %<>% str_replace_all(".*1.*", "Nodes from green subgraph")
set.seed(1)
left_annotation <- rowAnnotation(`Nodes` = node_types, col = list( `Nodes` = c(
  "Nodes from yellow subgraph"        = "gold", "Nodes from green subgraph"   = 'aquamarine3')))

right_annotation <-   
  rowAnnotation(gap = unit(12, "points"),Ce =  anno_barplot(bar_width = .5,width = unit(1, "cm"), border = T,
                                                           total.abs.centrality, gp = gpar(fill = 'azure4')))

ht_random <- Heatmap(corrM, name = "Correlation",  left_annotation =left_annotation, right_annotation=right_annotation,
        clustering_distance_columns = function(m)   dist(m, method = 'euclidean'),
        cluster_columns             = function(x) fastcluster::hclust(dist(x), "median"),
        clustering_distance_rows    = function(m)   dist(m, method = 'euclidean'),
        cluster_rows                = function(x) fastcluster::hclust(dist(x), "median"),
        row_km = 2, column_km = 2, border = TRUE,
        row_dend_width     = unit(.5, "cm"),
        column_dend_height = unit(.5,"cm"),
        show_column_names  = F,
        show_row_names     = F,
        row_names_gp = gpar(fontsize = 8),
        row_title_gp = gpar(fontsize = 8),
        column_title_gp = gpar(fontsize = 8),
        width = unit(6, "cm"), 
        height = unit(6, "cm"))

##########
#########

effects_on_sub0 <-  get_results("./results/0.xlsx", "fold_change_geometric","./results/0.txt", "./results/1.txt" )
effects_on_sub1 <-  get_results("./results/1.xlsx", "fold_change_geometric","./results/0.txt", "./results/1.txt" )

effects_on_sub0[[1]] %>% colnames() -> Ns
colnames(effects_on_sub0[[1]]) <- paste0(Ns, '_on_sub0')
effects_on_sub1[[1]] %>% colnames() -> Ns
colnames(effects_on_sub1[[1]]) <- paste0(Ns, '_on_sub1')

get_nodes_by_subgraph(effects_on_sub0[[1]], effects_on_sub0[[2]]) -> sub_0_nodes_by_subgraph
get_nodes_by_subgraph(effects_on_sub1[[1]], effects_on_sub1[[2]]) -> sub_1_nodes_by_subgraph

rbind(sub_0_nodes_by_subgraph$subgraph_0, sub_0_nodes_by_subgraph$subgraph_1) %>% as.data.frame() %>% 
  rename(subgraph_from_sub0 = subgraph) %>% rownames_to_column("node") -> effects_on_sub0_with_types

rbind(sub_1_nodes_by_subgraph$subgraph_0, sub_1_nodes_by_subgraph$subgraph_1) %>% as.data.frame()%>% 
  rename(subgraph_from_sub1 = subgraph)  %>% rownames_to_column("node")-> effects_on_sub1_with_types

first_join <- inner_join(effects_on_sub0_with_types, effects_on_sub1_with_types) 
first_join %>% column_to_rownames("node") %>% select_if(is.character) %>% .[[1]] -> node_types
first_join %>% column_to_rownames("node") %>% select_if(is.double)-> all_fc_centralities
rowSums(abs(all_fc_centralities)) -> total.abs.centrality



library(WGCNA)
all_fc_centralities %>% transposeBigData   %>% cor(method = "kendall") -> corrM
library(ComplexHeatmap)

node_types %<>% str_replace_all(".*0.*", "Nodes from yellow subgraph")
node_types %<>% str_replace_all(".*1.*", "Nodes from green subgraph")
set.seed(1)
left_annotation <- rowAnnotation(`Nodes` = node_types, col = list( `Nodes` = c(
  "Nodes from yellow subgraph"        = "gold", "Nodes from green subgraph"   = 'aquamarine3')))

right_annotation <-   
  rowAnnotation(gap = unit(12, "points"),Ce = anno_barplot(bar_width = .5,width = unit(1, "cm"), border = T,
                                                           total.abs.centrality, gp = gpar(fill = 'azure4')))


ht_freescale <- Heatmap(corrM, name = "Correlation",  left_annotation =left_annotation, right_annotation=right_annotation,
                     clustering_distance_columns = function(m)   dist(m, method = 'euclidean'),
                     cluster_columns             = function(x) fastcluster::hclust(dist(x), "median"),
                     clustering_distance_rows    = function(m)   dist(m, method = 'euclidean'),
                     cluster_rows                = function(x) fastcluster::hclust(dist(x), "median"),
                     row_km = 2, column_km = 2, border = TRUE,
                     row_dend_width     = unit(.5, "cm"),
                     column_dend_height = unit(.5,"cm"),
                     show_column_names  = F,
                     show_row_names     = F,
                     row_names_gp = gpar(fontsize = 8),
                     row_title_gp = gpar(fontsize = 8),
                     column_title_gp = gpar(fontsize = 8),
                     width = unit(6, "cm"), 
                     height = unit(6, "cm"))

library(grid)
grid.grabExpr(draw(ht_freescale)) -> ht_freescale
grid.grabExpr(draw(ht_random)) -> ht_random



bottom <- ggarrange(ht_freescale, ht_random,
          ncol = 2, nrow = 1,  labels = c('g','n'), widths = c(1,1), heights = c(1))



full_upper <- readRDS("./results/full_upper.rds")

full_panel <- ggarrange(full_upper, bottom, ncol = 1, nrow = 2, heights = c(1, .35))




ggsave(file="./results/centrality_panel.png",       plot=full_panel, width=11.5, height=15, dpi = 500)





