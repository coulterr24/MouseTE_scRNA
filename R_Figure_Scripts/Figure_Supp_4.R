
#### Figure Supp 4: Census of cell types of the mouse uterine tube ####

setwd("/Users/Coulter/Desktop/PhD Research/Seurat/Mouse Samples/Mouse_TE_Cancer")

#### Packages Load ####

library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(SoupX)
library(DoubletFinder)
library(data.table)
library(parallel)
library(tidyverse)
library(SoupX)
library(ggrepel)

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(viridisLite)
library(Polychrome)
library(circlize)
library(NatParksPalettes)

#### Proximal Datasets ####


Proximal <- readRDS( file = "../Proximal/20220907_Filtered_Cells_Final.rds" , refhook =  NULL)

Proximal_Named <- RenameIdents(Proximal, 
                               '0' = "Fibroblast 1", 
                               '1' = "Stem-like Epithelial", 
                               '2' = "Fibroblast 2",
                               '3' = "Fibroblast 3", 
                               '4' = "Immune", 
                               '5' = "Secretory Epithelial", 
                               '6' = "Endothelial",
                               '7' = "Ciliated Epithelial",
                               '8' = "Mesothelial", 
                               '9' = "Smooth Muscle")

Proximal_Named@active.ident <- factor(x = Proximal_Named@active.ident, 
                                      levels = c('Fibroblast 1',
                                                 'Fibroblast 2',
                                                 'Fibroblast 3',
                                                 'Smooth Muscle',
                                                 'Endothelial',
                                                 'Stem-like Epithelial',
                                                 'Secretory Epithelial',
                                                 'Ciliated Epithelial',
                                                 'Immune',
                                                 'Mesothelial'))

Proximal_Named <- SetIdent(Proximal_Named, value = Proximal_Named@active.ident)


Epi_Filter <- readRDS(file = "../Proximal/20220914_Proximal_Epi_Cells.rds" , refhook =  NULL)

Epi_Named <- RenameIdents(Epi_Filter, 
                          '0' = "Dbi+/Spdefhigh Secretory", 
                          '1' = "Bmpr1b+ Progenitor", 
                          '2' = "Wfdc2+ Secretory",
                          '3' = "Ciliated", 
                          '4' = "Sox17high Secretory", 
                          '5' = "Kcne3+ Secretory")

Epi_Named@active.ident <- factor(x = Epi_Named@active.ident, levels = c("Ciliated",
                                                                        "Dbi+/Spdefhigh Secretory",
                                                                        "Kcne3+ Secretory",
                                                                        "Sox17high Secretory",
                                                                        "Wfdc2+ Secretory",
                                                                        "Bmpr1b+ Progenitor"))


#### Figure Supp 4a: Proximal All Cell Types ####


Fibroblasts <- c('#FF9D00' , '#FFB653' , '#FFCB9A')   # Oranges
Muscle <- c('#E55451') # Reds
Endothelial <- c('#A0E6FF')  # Blues
Epi <-c('#6E3E6E','#CCCCFF','#DF73FF') # Purples
Immune <- c( '#5A5E6B' ) # Grey
Meso <- "#1F51FF" # Neon BLue

colors <- c(Fibroblasts, Muscle, Endothelial, Epi, Immune, Meso)


p1 <- DimPlot(
  Proximal_Named,
  reduction='umap',
  cols=colors,
  pt.size = 1.4,
  label.size = 4,
  label.color = "black",
  repel = TRUE,
  label=F) +
  NoLegend() +
  labs(x="UMAP_1",y="UMAP_2")


ggsave(filename = "FIGs3a_all_proximal_umap.pdf", plot = p1, width = 15, height = 12, dpi = 600)


#### Figure Supp 4b: Proximal Tile Mosaic  ####

library(treemap)

Proximal_Named <- RenameIdents(Proximal, 
                               '0' = "Fibroblast 1", 
                               '1' = "Stem-like Epithelial", 
                               '2' = "Fibroblast 2",
                               '3' = "Fibroblast 3", 
                               '4' = "Immune", 
                               '5' = "Secretory Epithelial", 
                               '6' = "Endothelial",
                               '7' = "Ciliated Epithelial",
                               '8' = "Mesothelial", 
                               '9' = "Smooth Muscle")

Proximal_Named@active.ident <- factor(x = Proximal_Named@active.ident, 
                                      levels = c('Fibroblast 1',
                                                 'Fibroblast 2',
                                                 'Fibroblast 3',
                                                 'Smooth Muscle',
                                                 'Endothelial',
                                                 'Stem-like Epithelial',
                                                 'Secretory Epithelial',
                                                 'Ciliated Epithelial',
                                                 'Immune',
                                                 'Mesothelial'))

Proximal_Named <- SetIdent(Proximal_Named, value = Proximal_Named@active.ident)



Fibroblasts <- c('#FF9D00' , '#FFB653' , '#FFCB9A')   # Oranges
Muscle <- c('#E55451') # Reds
Endothelial <- c('#A0E6FF')  # Blues
Epi <-c('#6E3E6E','#CCCCFF','#DF73FF') # Purples
Immune <- c( '#5A5E6B' ) # Grey
Meso <- "#1F51FF" # Neon BLue

colors <- c(Fibroblasts, Muscle, Endothelial, Epi, Immune, Meso)

prox_cell_types <- table(Idents(Proximal_Named), Proximal_Named$orig.ident)
prox_cell_type_df <- as.data.frame(prox_cell_types)


## Tile Mosaic ##

prox_treemap <- treemap(prox_cell_type_df, index = 'Var1', vSize= 'Freq', vColor = colors, palette = colors)


ggsave(filename = "20240612_all_prox_tile.pdf", plot = prox_cell_type_df, width = 8, height = 12, dpi = 600)


#### Figure Supp 4c: Epi Markers All Distal Cells ####



### Stacked Violin Plot Function ###

#https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat

modify_vlnplot <- function(obj, 
                           feature, 
                           pt.size = 0, 
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                           ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0, face = "bold.italic"), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 60, hjust=1, vjust=0.95), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  # plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
  plot_list<- purrr::map2(plot_list, c(5,5,8,5), function(x,y) x +                           
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

features<- c("Epcam", "Krt8" , "Ovgp1" , "Foxj1" )




Fibroblasts <- c('#FF9D00' , '#FFB653' , '#FFCB9A')   # Oranges
Muscle <- c('#E55451') # Reds
Endothelial <- c('#A0E6FF')  # Blues
Epi <-c('#6E3E6E','#CCCCFF','#DF73FF') # Purples
Immune <- c( '#5A5E6B' ) # Grey
Meso <- "#1F51FF" # Neon BLue

colors <- c(Fibroblasts, Muscle, Endothelial, Epi, Immune, Meso)

stack_vln <- StackedVlnPlot(obj = Proximal_Named, features = features, slot = "data",
                            pt.size = 0,
                            cols = colors)+ #9
  theme(plot.title = element_text(size = 32, face = "bold.italic"))+
  scale_x_discrete(limits = c('Fibroblast 1',
                              'Fibroblast 2',
                              'Fibroblast 3',
                              'Smooth Muscle',
                              'Endothelial',
                              'Stem-like Epithelial',
                              'Secretory Epithelial',
                              'Ciliated Epithelial',
                              'Immune',
                              'Mesothelial'))+
  theme(axis.text.x = element_text(size = 16, angle = 60))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.y.left = element_text(size = 16))

ggsave(filename = "20240612_All_Prox_stacked_vln.pdf", plot = stack_vln, width = 18, height = 12, dpi = 600)


#### Figure Supp 3d: Proximal Epi Cell Types ####


epi_umap <- DimPlot(object = Epi_Named,                # Seurat object  
                    reduction = 'umap',                 # Axes for the plot (UMAP, PCA, etc.) 
                    #group.by = "Patient",       # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
                    repel = TRUE,                       # Whether to repel the cluster labels
                    label = FALSE,                       # Whether to have cluster labels 
                    cols = c( "#35EFEF",
                              "#E95FE0",
                              "#B20224", 
                              "#F28D86", 
                              "#FB1111", 
                              "#FEB0DB"), 
                    
                    pt.size = 1.6,                      # Size of each dot is (0.1 is the smallest)
                    label.size = 2) +                   # Font size for labels    
  # You can add any ggplot2 1customizations here
  labs(title = 'Colored by Cluster')+        # Plot title
  NoLegend()
  
  
ggsave(filename = "FIGs3b_epi_proximal_umap.pdf", plot = epi_umap, width = 15, height = 12, dpi = 600)


#### Figure Supp 4e: Epi Markers All Distal Cells ####


P_Epi_Filter <- readRDS(file = "../Proximal/20220914_Proximal_Epi_Cells.rds" , refhook =  NULL)

P_Epi_Named <- RenameIdents(P_Epi_Filter, 
                            '0' = "Dbi+/Spdefhigh Secretory", 
                            '1' = "Bmpr1b+ Progenitor", 
                            '2' = "Wfdc2+ Secretory",
                            '3' = "Ciliated", 
                            '4' = "Sox17high Secretory", 
                            '5' = "Kcne3+ Secretory")

P_Epi_Named@active.ident <- factor(x = P_Epi_Named@active.ident, levels = c("Bmpr1b+ Progenitor",
                                                                            "Wfdc2+ Secretory",
                                                                            "Sox17high Secretory",
                                                                            "Kcne3+ Secretory",
                                                                            "Dbi+/Spdefhigh Secretory",
                                                                            "Ciliated"))



stack_vln <- StackedVlnPlot(obj = P_Epi_Named, features = features, slot = "data",
                            pt.size = 0,
                            cols = c( "#35EFEF",
                                      "#F28D86", 
                                      "#FB1111",
                                      "#FEB0DB",
                                      "#B20224",
                                      "#E95FE0"
                            ))+
  theme(plot.title = element_text(size = 32, face = "bold.italic"))+
  scale_x_discrete(limits = c("Bmpr1b+ Progenitor",
                              "Wfdc2+ Secretory",
                              "Sox17high Secretory",
                              "Kcne3+ Secretory",
                              "Dbi+/Spdefhigh Secretory",
                              "Ciliated"))+
  theme(axis.text.x = element_text(size = 16, angle = 60))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.y.left = element_text(size = 16))

ggsave(filename = "20240612_Prox_Epi_stacked_vln.pdf", plot = stack_vln, width = 18, height = 12, dpi = 600)


#### Figure Supp 4f: Proximal Epi Features####


named_features <- c("Krt8","Epcam", "Msln",
                    "Slc1a3","Sox9","Itga6", "Bmpr1b",
                    "Ovgp1","Sox17","Pax8", "Egr1",
                    "Wfdc2","Dbi","Gsto1","Fxyd4","Vim","Kcne3",
                    "Spdef","Lgals1","Upk1a", "Thrsp",
                    "Selenop", "Gstm2",
                    "Anpep", "Klf6",
                    "Id2",
                    "Ifit1",
                    "Prom1", "Ly6a", "Kctd8", "Adam8",
                    "Foxj1","Fam183b",
                    "Rgs22","Dnali1", "Mt1" , "Dynlrb2")



prox_dp <- DotPlot(object = Epi_Named,                    # Seurat object
                   assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                   features = named_features,                  # List of features (select one from above or create a new one)
                   # Colors to be used in the gradient
                   col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
                   col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
                   dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
                   dot.scale = 6,                        # Scale the size of the points
                   group.by = NULL,              # How the cells are going to be grouped
                   split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
                   scale = TRUE,                         # Whether the data is scaled
                   scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
                   scale.min = NA,                       # Set lower limit for scaling
                   scale.max = NA                        # Set upper limit for scaling
)+    labs(x = NULL,                              # x-axis label
           y = NULL)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 12 , face = "italic"))+
  theme(axis.text.y = element_text(size = 12))+
  theme(legend.title = element_text(size = 12))+ 
  scale_y_discrete(limits = c("Ciliated",
                              "Dbi+/Spdefhigh Secretory",
                              "Kcne3+ Secretory",
                              "Sox17high Secretory",
                              "Wfdc2+ Secretory",
                              "Bmpr1b+ Progenitor"
  ))

ggsave(filename = "FIGs3c_epi_proximal_dotplot.pdf", plot = prox_dp, width = 18, height = 10, dpi = 600)



