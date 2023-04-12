library(tidyverse)
library(ggplot2)
library(ggpubr)
get_results <- function(xlsx_path, my_sheet, nodes_of_subgraph_0_path , nodes_of_subgraph_1_path){
  
  df1 <- readxl::read_xlsx(xlsx_path, sheet = my_sheet) 
  colnames(df1) -> my_cols
  
  df2                        <- df1 %>% column_to_rownames(my_cols[1])
  subgraph_0_nodes           <- read.delim(nodes_of_subgraph_0_path, header = F)[[1]] %>% c
  subgraph_1_nodes           <-  read.delim(nodes_of_subgraph_1_path, header = F)[[1]] %>% c
  return(list(df2, subgraph_0_nodes, subgraph_1_nodes))}


get_nodes_by_subgraph <- function(centality_variations, subgraph_0_nodes){
  centality_variations %>% mutate(subgraph = if_else(row.names(.) %in% subgraph_0_nodes, "subgraph_0_node", "subgraph_1_node" ) ) -> df1
  df1 %>% filter(subgraph == "subgraph_0_node") -> subgraph_0
  df1 %>%  filter(subgraph == "subgraph_1_node") -> subgraph_1
  list(subgraph_0, subgraph_1) %>% purrr::set_names(c('subgraph_0', 'subgraph_1')) -> nodes_by_subgraph
  return(nodes_by_subgraph)
}

plot_density <- function(df, titulo='sin título', fill = "deepskyblue4" , xlab = "Centrality", ylab = "# nodes", xlim = c(-0.04, 0.04)){cut = 0.005
df %>% gather() ->df_gathered
df_gathered['value'] %>% filter(abs(value) > cut) %>%
  ggplot( aes(x=value)) +
  geom_density(fill=fill, color="azure4", alpha=.8)+
  theme(legend.position="top" , plot.title = element_text(size=8.5)) +
  ylab(ylab) + ggtitle(titulo) +  xlim(xlim) +
  xlab(xlab) + geom_vline(xintercept=0, linetype="dashed", color = "black",  size=1) -> Density_plot_of_centralities
return(Density_plot_of_centralities)
}

generate_subplots <- function(nodes_vs_subsystems, color = 'brown4'){
  
  get_a_subplot <- function(nodes_vs_subsystems,a_effect,xlab,color ){
    pluck(nodes_vs_subsystems, a_effect)  -> df
    df %>% min -> xmin
    df %>% max -> xmax 
    df %>% plot_density(titulo =a_effect, fill = color, xlab = xlab, xlim = c(xmin*1.2, xmax*1.2)) }
  
  get_a_subplot(nodes_vs_subsystems, a_effect = "Yellow subgraph upon yellow removals",color = color, xlab = "Centrality fold-change") -> A
  get_a_subplot(nodes_vs_subsystems, a_effect = 'effect_on_sub0_from_sub1_nodes',color = color, xlab = "Centrality fold-change") -> B
  get_a_subplot(nodes_vs_subsystems, a_effect = 'effect_on_sub1_from_sub0_nodes',color = color, xlab = "Centrality fold-change") -> C
  get_a_subplot(nodes_vs_subsystems, a_effect = 'effect_on_sub1_from_sub1_nodes',color = color, xlab = "Centrality fold-change") -> D
  return(list(A,B,C,D) %>% 
           purrr::set_names(c('A','B','C','D')))}
