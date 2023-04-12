library(ggplot2)
library(ggpubr)
img1 <- png::readPNG("./Supplementary/results/model_rxns.png")
img2 <-  png::readPNG("./Supplementary/results/phpp_sampleo.png")
img3 <-  png::readPNG("./Supplementary/results/phpp_fba_solution.png")
img4 <-  png::readPNG("./Supplementary/results/Bipartite_small_metabolism.png")

im_a <- ggplot() + background_image(img1) # tried to insert the image as background.. there must be a better way
im_b <- ggplot() + background_image(img2)
im_c <- ggplot() + background_image(img3) # tried to insert the image as background.. there must be a better way
im_d <- ggplot() + background_image(img4) 



upper <- ggarrange(im_a,im_b, im_c, ncol = 3, nrow = 1, widths = c(.8,.7,.7) ,
                   labels = c("a", "b","c"),vjust = .9, hjust = 0.1, align = "hv" )

panel <- ggarrange(upper,im_d ,ncol = 1, nrow = 2 , heights = c(0.4,1) , labels = c("","d"),vjust = 1, hjust = 0.1)


ggsave(file="./Supplementary/results/toy_metabolism_phpp_fba.png",       plot=panel, width=7, height=8, dpi = 500)
