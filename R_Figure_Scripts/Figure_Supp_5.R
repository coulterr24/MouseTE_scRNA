#### Figure Supp 5: Distal and proximal epithelial cell correlation ####

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



#### Proximal vs Distal Cluster Correlation ####


Distal_Epi_Filter <- readRDS(file = "../Distal/20220817_Distal_Epi_Cells.rds" , refhook =  NULL)

Distal_Epi_Named <- RenameIdents(Distal_Epi_Filter, 
                                 '0' = "Spdef+ Secretory", 
                                 '1' = "Slc1a3+ Stem/Progenitor", 
                                 '2' = "Cebpdhigh/Foxj1- Progenitor",
                                 '3' = "Ciliated 1", 
                                 '4' = "Ciliated 2", 
                                 '5' = "Pax8low/Prom1+ Cilia-forming", 
                                 '6' = "Fibroblast-like",
                                 '7' = "Slc1a3med/Sox9+ Cilia-forming",
                                 '8' = "Selenop+/Gstm2high Secretory")

Distal_Epi_Named@active.ident <- factor(x = Distal_Epi_Named@active.ident, levels = c( c("Slc1a3+ Stem/Progenitor",
                                                                                         "Cebpdhigh/Foxj1- Progenitor",
                                                                                         "Slc1a3med/Sox9+ Cilia-forming",
                                                                                         "Pax8low/Prom1+ Cilia-forming", 
                                                                                         "Fibroblast-like",
                                                                                         "Spdef+ Secretory",
                                                                                         "Selenop+/Gstm2high Secretory",
                                                                                         "Ciliated 1",
                                                                                         "Ciliated 2")))



prox_avg_exp <- AverageExpression(Epi_Named)$RNA
distal_avg_exp <- AverageExpression(Distal_Epi_Named)$RNA


cor.exp <- as.data.frame(cor(x = prox_avg_exp , y = distal_avg_exp))

cor.exp$x <- rownames(cor.exp)

cor.df <- tidyr::gather(data = cor.exp, y, correlation, c("Slc1a3+ Stem/Progenitor",
                                                          "Cebpdhigh/Foxj1- Progenitor",
                                                          "Slc1a3med/Sox9+ Cilia-forming",
                                                          "Pax8low/Prom1+ Cilia-forming", 
                                                          "Fibroblast-like",
                                                          "Spdef+ Secretory",
                                                          "Selenop+/Gstm2high Secretory",
                                                          "Ciliated 1",
                                                          "Ciliated 2"))


distal_cells <- c("Slc1a3+ Stem/Progenitor",
                  "Cebpdhigh/Foxj1- Progenitor",
                  "Slc1a3med/Sox9+ Cilia-forming",
                  "Pax8low/Prom1+ Cilia-forming", 
                  "Fibroblast-like",
                  "Spdef+ Secretory",
                  "Selenop+/Gstm2high Secretory",
                  "Ciliated 1",
                  "Ciliated 2")

prox_cells <- c("Bmpr1b+ Progenitor",
                "Ciliated",
                "Dbi+/Spdefhigh Secretory",
                "Wfdc2+ Secretory",
                "Sox17high Secretory",
                "Kcne3+ Secretory")


corr_matrix <- ggplot(cor.df, aes(x, y, fill = correlation)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_viridis_c(values = c(0,1),option="rocket", begin=.4,end=0.99, direction = -1,)+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 12, face = "bold"),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 12, face = "bold.italic"))+
  theme(plot.title = element_blank())+
  scale_y_discrete(limits = c("Ciliated 2",
                              "Ciliated 1",
                              "Selenop+/Gstm2high Secretory",
                              "Spdef+ Secretory",
                              "Fibroblast-like",
                              "Pax8low/Prom1+ Cilia-forming", 
                              "Slc1a3med/Sox9+ Cilia-forming",
                              "Cebpdhigh/Foxj1- Progenitor",
                              "Slc1a3+ Stem/Progenitor"))+
  scale_x_discrete(limits = c("Bmpr1b+ Progenitor",
                              "Wfdc2+ Secretory",
                              "Sox17high Secretory",
                              "Kcne3+ Secretory",
                              "Dbi+/Spdefhigh Secretory",
                              "Ciliated"))+
  geom_text(aes(x, y, label = round(correlation, digits = 2)), color = "black", size = 4)


ggsave(filename = "FIGs3d_epi_cluster_corr.pdf", plot = corr_matrix, width = 18, height = 10, dpi = 600)




