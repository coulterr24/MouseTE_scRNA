
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


#### Load Distal and Proximal Epithelial Datasets ####

Distal_Epi_Filter <- readRDS(file = "../Distal/20220817_Distal_Epi_Cells.rds" , refhook =  NULL)

Distal_Epi_Filter <- RenameIdents(Distal_Epi_Filter, 
                          '0' = "Spdef+ Secretory", 
                          '1' = "Slc1a3+ Stem/Progenitor", 
                          '2' = "Cebpdhigh/Foxj1- Progenitor",
                          '3' = "Ciliated 1", 
                          '4' = "Ciliated 2", 
                          '5' = "Pax8low/Prom1+ Cilia-forming", 
                          '6' = "Fibroblast-like",
                          '7' = "Slc1a3med/Sox9+ Cilia-forming",
                          '8' = "Selenop+/Gstm2high Secretory")

Distal_Epi_Filter@active.ident <- factor(x = Distal_Epi_Filter@active.ident, levels = c( c("Slc1a3+ Stem/Progenitor",
                                                                           "Cebpdhigh/Foxj1- Progenitor",
                                                                           "Slc1a3med/Sox9+ Cilia-forming",
                                                                           "Pax8low/Prom1+ Cilia-forming", 
                                                                           "Fibroblast-like",
                                                                           "Spdef+ Secretory",
                                                                           "Selenop+/Gstm2high Secretory",
                                                                           "Ciliated 1",
                                                                           "Ciliated 2")))





Prox_Epi_Filter <- readRDS(file = "../Proximal/20220914_Proximal_Epi_Cells.rds" , refhook =  NULL)

Prox_Epi_Filter <- RenameIdents(Prox_Epi_Filter, 
                          '0' = "Dbi+/Spdefhigh Secretory", 
                          '1' = "Bmpr1b+ Progenitor", 
                          '2' = "Wfdc2+ Secretory",
                          '3' = "Ciliated", 
                          '4' = "Sox17high Secretory", 
                          '5' = "Kcne3+ Secretory")

Prox_Epi_Filter@active.ident <- factor(x = Prox_Epi_Filter@active.ident, levels = c("Ciliated",
                                                                        "Dbi+/Spdefhigh Secretory",
                                                                        "Kcne3+ Secretory",
                                                                        "Sox17high Secretory",
                                                                        "Wfdc2+ Secretory",
                                                                        "Bmpr1b+ Progenitor"))



#### Retrieving Pax8+ Distal Epithelial Cell Percentages ####


## PAX8+ Cells ##


plot <- DotPlot(object = Distal_Epi_Named, features = c("Pax8"))

plot_data <- plot$data %>% 
  select(pct.exp, id)

## Cell Number per cluster ##

cells <- Distal_Epi_Named@active.ident %>% as.data.table
cluster_cell_number <- cells[, .N, by = cells[,]]

cluster_cell_number <- cluster_cell_number %>% rename(. = 'id')

## Combined Table ##

dis_combined <- merge(plot_data,cluster_cell_number)

write.csv(dis_combined, file = "TBLs1_Distal_Pax8_Percent.csv")



#### Retrieving Pax8+ Proximal Epithelial Cell Percentages ####


plot <- DotPlot(object = Prox_Epi_Named, features = c("Pax8"))

plot_data <- plot$data %>% 
  select(pct.exp, id)

## Cell Number per cluster ##

cells <- Prox_Epi_Named@active.ident %>% as.data.table
cluster_cell_number <- cells[, .N, by = cells[,]]

cluster_cell_number <- cluster_cell_number %>% rename(. = 'id')

## Combined Table ##

prox_combined <- merge(plot_data,cluster_cell_number)

write.csv(prox_combined, file = "TBLs1_Prox_Pax8_Percent.csv")
