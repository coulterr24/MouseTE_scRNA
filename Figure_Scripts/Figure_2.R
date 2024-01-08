
#### Figure 2: Characterization of distal epithelial cell states ####

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


#### Load Distal Epithelial Dataset ####

Epi_Filter <- readRDS(file = "../Distal/20220817_Distal_Epi_Cells.rds" , refhook =  NULL)

Epi_Named <- RenameIdents(Epi_Filter, 
                          '0' = "Spdef+ Secretory", 
                          '1' = "Slc1a3+ Stem/Progenitor", 
                          '2' = "Cebpdhigh/Foxj1- Progenitor",
                          '3' = "Ciliated 1", 
                          '4' = "Ciliated 2", 
                          '5' = "Pax8low/Prom1+ Cilia-forming", 
                          '6' = "Fibroblast-like",
                          '7' = "Slc1a3med/Sox9+ Cilia-forming",
                          '8' = "Selenop+/Gstm2high Secretory")

Epi_Named@active.ident <- factor(x = Epi_Named@active.ident, levels = c( c("Slc1a3+ Stem/Progenitor",
                                                                           "Cebpdhigh/Foxj1- Progenitor",
                                                                           "Slc1a3med/Sox9+ Cilia-forming",
                                                                           "Pax8low/Prom1+ Cilia-forming", 
                                                                           "Fibroblast-like",
                                                                           "Spdef+ Secretory",
                                                                           "Selenop+/Gstm2high Secretory",
                                                                           "Ciliated 1",
                                                                           "Ciliated 2")))



#### Figure 2a: Epithelial Cells of the Distal Uterine Tube ####

epi_umap <- DimPlot(object = Epi_Named,                # Seurat object  
                    reduction = 'umap',                 # Axes for the plot (UMAP, PCA, etc.) 
                    repel = TRUE,                       # Whether to repel the cluster labels
                    label = FALSE,                       # Whether to have cluster labels 
                    cols = c("#35EFEF", #1
                             "#00A1C6", #2
                             "#2188F7", #3
                             "#EA68E1", #4
                             "#59D1AF", #5
                             "#B20224", #6
                             "#F28D86", #7
                             "#A374B5", #8
                             "#9000C6"), #9
                    
                    pt.size = 0.6,                      # Size of each dot is (0.1 is the smallest)
                    label.size = 0.5) +                   # Font size for labels    
  # You can add any ggplot2 1customizations here
  labs(title = 'Colored by Cluster')+        # Plot title
  NoLegend() +
  labs(x="UMAP_1",y="UMAP_2")

ggsave(filename = "Fig2a_epi_umap.pdf", plot = epi_umap, width = 15, height = 12, dpi = 600)


#### Figure 2b: Distal Uterine Tube Features for Epithelial Cell State Identification ####

distal_features <- c("Krt8","Epcam",
                     "Slc1a3","Cd44","Sox9",
                     "Ovgp1","Sox17","Pax8", "Egr1",
                     "Itga6", "Bmpr1b",
                     "Rhoj", "Klf6","Msln","Cebpd",
                     "Dpp6", "Sec14l3", "Fam161a",
                     "Prom1", "Ly6a", "Kctd8", "Adam8",
                     "Dcn", "Col1a1", "Col1a2", "Timp3", "Pdgfra","Lgals1",
                     "Upk1a", "Thrsp","Spdef",
                     "Selenop", "Gstm2",
                     "Foxj1","Fam183b",
                     "Rgs22","Dnali1", "Mt1" , "Dynlrb2")


epi_dp <- DotPlot(object = Epi_Named,                    # Seurat object
                  assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                  features = distal_features,                  # List of features (select one from above or create a new one)
                  # Colors to be used in the gradient
                  col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
                  col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
                  dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
                  dot.scale = 4,                        # Scale the size of the points
                  group.by = NULL,              # How the cells are going to be grouped
                  split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
                  scale = TRUE,                         # Whether the data is scaled
                  scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
                  scale.min = NA,                       # Set lower limit for scaling
                  scale.max = NA )+                       # Set upper limit for scaling
  labs(x = NULL,                              # x-axis label
       y = NULL)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  #theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+
  theme(axis.text.x = element_text(size = 8 , face = "italic"))+
  theme(axis.text.y = element_text(size = 9))+
  theme(legend.title = element_text(size = 9))+
  theme(legend.text = element_text(size = 8))+ 
  scale_y_discrete(limits = c("Ciliated 2",
                              "Ciliated 1",
                              "Selenop+/Gstm2high Secretory",
                              "Spdef+ Secretory",
                              "Fibroblast-like",
                              "Pax8low/Prom1+ Cilia-forming", 
                              "Slc1a3med/Sox9+ Cilia-forming",
                              "Cebpdhigh/Foxj1- Progenitor",
                              "Slc1a3+ Stem/Progenitor"))

ggsave(filename = "Fig2b_epi_dot_plot.pdf", plot = epi_dp, width = 8.3, height = 4.0625, dpi = 600)













