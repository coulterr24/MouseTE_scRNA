
#### Figure Supp 3: Census of cell types of the mouse uterine tube ####

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


#### Figure Supp 3a: Proximal All Cell Types ####


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


#### Figure Supp 3b: Proximal Epi Cell Types ####


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


#### Figure Supp 3c: Proximal Epi Features####


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


#### Figure Supp 3d: Proximal vs Distal Cluster Correlation ####


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




