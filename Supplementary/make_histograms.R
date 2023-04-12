setwd("~/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism/Supplementary")
source('R_functions.R')
#get data


effects_on_sub0 <-  get_results("./results/proj_two_core_diRandom_0.xlsx", "fold_change_geometric","./results/proj_two_core_diRandom_0.txt", "./results/proj_two_core_diRandom_1.txt" )
effects_on_sub1 <-  get_results("./results/proj_two_core_diRandom_1.xlsx", "fold_change_geometric","./results/proj_two_core_diRandom_0.txt", "./results/proj_two_core_diRandom_1.txt" )
#Manipulations
get_nodes_by_subgraph(effects_on_sub0[[1]], effects_on_sub0[[2]]) -> sub_0_nodes_by_subgraph
get_nodes_by_subgraph(effects_on_sub1[[1]], effects_on_sub1[[2]]) -> sub_1_nodes_by_subgraph

sub_0_nodes_by_subgraph$subgraph_0  %>% select(-"subgraph") -> `Yellow subgraph upon yellow removals`
sub_0_nodes_by_subgraph$subgraph_1  %>% select(-"subgraph")-> effect_on_sub0_from_sub1_nodes
sub_1_nodes_by_subgraph$subgraph_0  %>% select(-"subgraph")-> effect_on_sub1_from_sub0_nodes
sub_1_nodes_by_subgraph$subgraph_1  %>% select(-"subgraph")-> effect_on_sub1_from_sub1_nodes

subgraphs_vs_nodes <- list(`Yellow subgraph upon yellow removals`, effect_on_sub0_from_sub1_nodes, effect_on_sub1_from_sub0_nodes, effect_on_sub1_from_sub1_nodes) %>% 
                      purrr::set_names('Yellow subgraph upon yellow removals', 'Yellow subgraph upon green removals', 
                                       'Green subgraph upon yellow removals', 'Green subgraph upon green removals') 


generate_subplots <- function(nodes_vs_subsystems, color = 'brown4'){
  get_a_subplot <- function(nodes_vs_subsystems,a_effect,xlab,color ){
    pluck(nodes_vs_subsystems, a_effect)  -> df
    df %>% min -> xmin
    df %>% max -> xmax 
    df %>% plot_density(titulo =a_effect, fill = color, xlab = xlab, xlim = c(xmin*1.2, xmax*1.2)) }
  
  get_a_subplot(nodes_vs_subsystems, a_effect = "Yellow subgraph upon yellow removals",color = color, xlab = "Centrality fold-change") -> A
  get_a_subplot(nodes_vs_subsystems, a_effect = 'Yellow subgraph upon green removals',color = color, xlab = "Centrality fold-change") -> B
  get_a_subplot(nodes_vs_subsystems, a_effect = 'Green subgraph upon yellow removals',color = color, xlab = "Centrality fold-change") -> C
  get_a_subplot(nodes_vs_subsystems, a_effect = 'Green subgraph upon green removals',color = color, xlab = "Centrality fold-change") -> D
  return(list(A,B,C,D) %>% 
           purrr::set_names(c('A','B','C','D')))}



generate_subplots(subgraphs_vs_nodes, color = 'purple4') -> sublots_lists

########################################################3

panel_random     <- ggarrange(sublots_lists[[1]],sublots_lists[[2]], sublots_lists[[3]], 
                              sublots_lists[[4]], ncol = 1, nrow = 4,  
                              labels = c('j','k','l','m'), widths = c(1,1), heights = c(1,1))

img1 <- png::readPNG("./results/old_two_core_diRandom.png")
img2 <-  png::readPNG("./results/old_Proj_two_core_diRandom.png")

im_A <- ggplot() + background_image(img1) # tried to insert the image as background.. there must be a better way
im_B <- ggplot() + background_image(img2) 

top_panel <- ggarrange(im_A, im_B, labels = c("h", "i") )

upper_panel_right <-  ggarrange(top_panel, panel_random, ncol = 1, nrow = 2, heights = c(.4,1))

####
########

effects_on_sub0 <-  get_results("./results/0.xlsx", "fold_change_geometric","./results/0.txt", "./results/1.txt" )
effects_on_sub1 <-  get_results("./results/1.xlsx", "fold_change_geometric","./results/0.txt", "./results/1.txt" )
#Manipulations
get_nodes_by_subgraph(effects_on_sub0[[1]], effects_on_sub0[[2]]) -> sub_0_nodes_by_subgraph
get_nodes_by_subgraph(effects_on_sub1[[1]], effects_on_sub1[[2]]) -> sub_1_nodes_by_subgraph

sub_0_nodes_by_subgraph$subgraph_0  %>% select(-"subgraph") -> effect_on_sub0_from_sub0_nodes
sub_0_nodes_by_subgraph$subgraph_1  %>% select(-"subgraph")-> effect_on_sub0_from_sub1_nodes
sub_1_nodes_by_subgraph$subgraph_0  %>% select(-"subgraph")-> effect_on_sub1_from_sub0_nodes
sub_1_nodes_by_subgraph$subgraph_1  %>% select(-"subgraph")-> effect_on_sub1_from_sub1_nodes
subgraphs_vs_nodes <- list(effect_on_sub0_from_sub0_nodes, effect_on_sub0_from_sub1_nodes, effect_on_sub1_from_sub0_nodes, effect_on_sub1_from_sub1_nodes) %>% 
  purrr::set_names('Yellow subgraph upon yellow removals', 'Yellow subgraph upon green removals', 
                   'Green subgraph upon yellow removals', 'Green subgraph upon green removals') 

generate_subplots(subgraphs_vs_nodes, color = 'brown4') -> sublots_lists

panel_free_scale <- ggarrange(sublots_lists[[1]],sublots_lists[[2]], sublots_lists[[3]], sublots_lists[[4]], 
                              ncol = 1, nrow = 4,  labels = c('c','d','e','f'), widths = c(1,1), heights = c(1,1))


#############
img1 <- png::readPNG("./results/old_two_core_diFreeScale.png")
img2 <-  png::readPNG("./results/old_Proj_two_core_diFreeScale.png")

im_A <- ggplot() + background_image(img1) # tried to insert the image as background.. there must be a better way
im_B <- ggplot() + background_image(img2) 

top_panel <- ggarrange(im_A, im_B, labels = c("a", "b") )

upper_panel_left <-  ggarrange(top_panel, panel_free_scale, ncol = 1, nrow = 2, heights = c(.4,1))

full_upper <- ggarrange(upper_panel_left, upper_panel_right, ncol = 2, nrow = 1)



saveRDS(full_upper, "./results/full_upper.rds")




