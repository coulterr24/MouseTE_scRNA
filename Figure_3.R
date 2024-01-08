
#### Figure 3: Feature expression over the PHATE embedding ####

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

library(monocle3)

#### Load Distal Epithelial Pseudotime Dataset ####

Distal_PHATE <- readRDS(file = "../Distal/20220826_Distal_Epi_PHATE.rds" , refhook = NULL)

Beeg_PHATE <- Distal_PHATE

Beeg_PHATE@reductions[["phate"]]@cell.embeddings <- Distal_PHATE@reductions[["phate"]]@cell.embeddings*100

cds <- readRDS(file = "20221101_Distal_Epi_PHATE_Monocle3.rds" , refhook = NULL)


#### Figure 3a: PHATE embedding for differentiation trajectory of distal epithelial cells ####


phate_dif <- DimPlot(Beeg_PHATE , 
                    reduction = "phate", 
                    cols = c("#B20224", 
                             "#35EFEF", 
                             "#00A1C6", 
                             "#A374B5", 
                             "#9000C6", 
                             "#EA68E1", 
                             "#59D1AF", 
                             "#2188F7", 
                             "#F28D86"),
                    pt.size = 0.7,
                    shuffle = TRUE,
                    seed = 0,
                    label = FALSE)+  
              labs(title = 'Colored by Cluster')+        # Plot title
              NoLegend() +
              labs(x="UMAP_1",y="UMAP_2")

ggsave(filename = "Fig3a_epi_phate.pdf", plot = phate_dif, width = 15, height = 12, dpi = 600)


#### Figure 3b: PHATE and Monocle3 differentiation trajectory path ####


pseudtotime <- plot_cells(cds, 
                          color_cells_by = "pseudotime",
                          label_cell_groups=FALSE,
                          label_leaves=FALSE,
                          label_branch_points=FALSE,
                          graph_label_size=0,
                          cell_size = .01,
                          cell_stroke = 1)+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())+
  theme(legend.text = element_text(size = 12))

ggsave(filename = "Fig3b_epi_pseudtotime.pdf", plot = pseudtotime, width = 18, height = 12, dpi = 600)



#### Figure 3c: Slc1a3 PHATE Feature Plot ####


Slc1a3_PHATE <- FeaturePlot(Beeg_PHATE, features = c("Slc1a3"), reduction = "phate", pt.size = 0.5)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())+
  theme(legend.text = element_text(size = 12))

ggsave(filename = "Fig3c_SLC1A3_PHATE.pdf", plot = Slc1a3_PHATE, width = 18, height = 9, dpi = 600)


#### Figure 3d: Pax8 PHATE Feature Plot ####

Pax8_PHATE <- FeaturePlot(Beeg_PHATE, features = c("Pax8"), reduction = "phate", pt.size = 0.5)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())+
  theme(legend.text = element_text(size = 12))

ggsave(filename = "Fig3d_PAX8_PHATE.pdf", plot = Pax8_PHATE, width = 18, height = 9, dpi = 600)

