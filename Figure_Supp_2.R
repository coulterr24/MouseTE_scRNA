
#### Figure Supp 2: Characterization of distal epithelial cell states ####

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


#### Unprocessed and Processed Distal Dataset ####

All.merged <- readRDS(file = "../Distal/20220810_Mouse_Distal.rds", refhook = NULL) # Prior to Quality Control

Distal <- readRDS(file = "../Distal/20220810_Filtered_Cells.rds" , refhook =  NULL) # After to Quality Control

Distal_Named <- RenameIdents(Distal, 
                             '0' = "Fibroblast 1", 
                             '1' = "Fibroblast 2", 
                             '2' = "Secretory Epithelial",
                             '3' = "Smooth Muscle", 
                             '4' = "Ciliated Epithelial 1", 
                             '5' = "Fibroblast 3", 
                             '6' = "Stem-like Epithelial 1",
                             '7' = "Stem-like Epithelial 2",
                             '8' = "Ciliated Epithelial 2", 
                             '9' = "Blood Endothelial", 
                             '10' = "Pericyte", 
                             '11' = "Intermediate Epithelial", 
                             '12' = "T/NK Cell", 
                             '13' = "Epithelial/Fibroblast", 
                             '14' = "Macrophage", 
                             '15' = "Erythrocyte", 
                             '16' = "Luteal",
                             '17' = "Mesothelial",
                             '18' = "Lymph Endothelial/Epithelial") # Remove cluster due few data points and suspected doublet

Distal_Named@active.ident <- factor(x = Distal_Named@active.ident, 
                                    levels = c('Fibroblast 1',
                                               'Fibroblast 2',
                                               'Fibroblast 3',
                                               'Smooth Muscle',
                                               'Pericyte',
                                               'Blood Endothelial',
                                               'Lymph Endothelial/Epithelial',
                                               'Epithelial/Fibroblast',
                                               'Stem-like Epithelial 1',
                                               'Stem-like Epithelial 2',
                                               'Intermediate Epithelial',
                                               'Secretory Epithelial',
                                               'Ciliated Epithelial 1',
                                               'Ciliated Epithelial 2',
                                               'T/NK Cell',
                                               'Macrophage',
                                               'Erythrocyte',
                                               'Mesothelial',
                                               'Luteal'))


Distal_Named <- subset(Distal_Named, 
                       idents = c('Fibroblast 1',
                                  'Fibroblast 2',
                                  'Fibroblast 3',
                                  'Smooth Muscle',
                                  'Pericyte',
                                  'Blood Endothelial',
                                  'Epithelial/Fibroblast',
                                  'Stem-like Epithelial 1',
                                  'Stem-like Epithelial 2',
                                  'Intermediate Epithelial',
                                  'Secretory Epithelial',
                                  'Ciliated Epithelial 1',
                                  'Ciliated Epithelial 2',
                                  'T/NK Cell',
                                  'Macrophage',
                                  'Erythrocyte',
                                  'Mesothelial',
                                  'Luteal'))

Distal_Named <- SetIdent(Distal_Named, value = Distal_Named@active.ident)


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

#### Figure Supp 2a: Unfilitered % MT Genes ####

unfiltered_MT <- VlnPlot(All.merged, features = c("percent.mt"), group.by = 'Sample', pt.size = 0,
                         cols = natparks.pals(name="Arches2",n=3))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+     # Change X-Axis Text Size
  theme(axis.text.y = element_text(size = 16))+     # Change Y-Axis Text Size
  theme(axis.title.y = element_text(size = 18))+    # Change Y-Axis Title Text Size
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank())

ggsave(filename = "FIGs2a_unfiltered_MT.pdf", plot = unfiltered_MT, width = 12, height = 9, dpi = 600)

#### Figure Supp 2b: Unfilitered nFeature RNA #### 

unfiltered_nFeature <- VlnPlot(All.merged, features = c("nFeature_RNA"), group.by = 'Sample', pt.size = 0,
                               cols = natparks.pals(name="Arches2",n=3))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 18))+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank()) # Change object to visualize other samples

ggsave(filename = "FIGs2b_unfiltered_nFeature.pdf", plot = unfiltered_nFeature, width = 12, height = 9, dpi = 600)

#### Figure Supp 2c: Unfilitered nCount RNA #### 

unfiltered_nCount <- VlnPlot(All.merged, features = c("nCount_RNA"), group.by = 'Sample', pt.size = 0,
                             cols = natparks.pals(name="Arches2",n=3))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 18))+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank()) # Change object to visualize other samples

ggsave(filename = "FIGs2c_unfiltered_nCount.pdf", plot = unfiltered_nCount, width = 12, height = 9, dpi = 600)

#### Figure Supp 2d: Doublets All Cells #### 

All_Doublet <- DimPlot(object = Distal, 
                       reduction = 'umap', 
                       group.by = "Doublet",
                       cols = c( "#ffb6c1", "#380b11"),
                       repel = TRUE, 
                       label = F, 
                       pt.size = 1.2, 
                       order = c("Doublet","Singlet"),
                       label.size = 5) +
  labs(x="UMAP_1",y="UMAP_2")

ggsave(filename = "FIGs2d1_All_Doublet_umap.pdf", plot = All_Doublet, width = 22, height = 17, dpi = 600)

## Stacked Bar Doublets ##

table <- table(Distal_Named@active.ident ,
               Distal_Named@meta.data$Doublet)    # Create a table of counts

df <- data.frame(table) 



doublet <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
                      aes(x = Var1,              # Variable to plot on the x-axis
                          y = Freq,              # Variable to plot on the y-axis
                          fill = factor(Var2,    # Variable to fill the bars
                                        levels = c("Doublet","Singlet")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Doublet", limits = c("Doublet","Singlet"),
                    values = c('#8B0000','#808080')) +

  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))+       # Text color and horizontal adjustment on y-axis
  scale_x_discrete(limits = (c('Intermediate Epithelial',
                               'Epithelial/Fibroblast',
                               'Stem-like Epithelial 1',
                               'Ciliated Epithelial 1',
                               'Erythrocyte',
                               'Smooth Muscle',
                               'Stem-like Epithelial 2',
                               'Mesothelial',
                               'Blood Endothelial',
                               'Pericyte',
                               'Fibroblast 2',
                               'Secretory Epithelial',
                               'Fibroblast 1',
                               'Ciliated Epithelial 2',
                               'Fibroblast 3',
                               'T/NK Cell',
                               'Macrophage',
                               'Luteal')))+
  coord_flip()



ggsave(filename = "FIGs2d2_doublet_quant.pdf", plot = doublet, width = 10, height = 16, dpi = 600)


#### Figure Supp 2e: Doublets Epithelial Cells #### 


Epi_Doublet <- DimPlot(object = Epi_Named, 
                       reduction = 'umap', 
                       group.by = "Doublet",
                       cols = c( "#ffb6c1", "#380b11"),
                       repel = TRUE, 
                       label = F, 
                       pt.size = 1.2, 
                       order = c("Doublet","Singlet"),
                       label.size = 5) +
  labs(x="UMAP_1",y="UMAP_2")

ggsave(filename = "FIGs2e1_Epi_Doublet_umap.pdf", plot = All_Doublet, width = 22, height = 17, dpi = 600)


## Stacked Bar Doublets ##

table <- table(Epi_Named@active.ident ,
               Epi_Named@meta.data$Doublet)    # Create a table of counts

df <- data.frame(table) 




epi_doublet <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
                      aes(x = Var1,              # Variable to plot on the x-axis
                          y = Freq,              # Variable to plot on the y-axis
                          fill = factor(Var2,    # Variable to fill the bars
                                        levels = c("Doublet","Singlet")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Doublet", limits = c("Doublet","Singlet"),
                    values = c('#8B0000','#808080')) +

  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))+       # Text color and horizontal adjustment on y-axis
  scale_x_discrete(limits = (c("Slc1a3med/Sox9+ Cilia-forming",
                               "Fibroblast-like",
                               "Pax8low/Prom1+ Cilia-forming",
                               "Slc1a3+ Stem/Progenitor",
                               "Cebpdhigh/Foxj1- Progenitor",
                               "Ciliated 1",
                               "Ciliated 2", 
                               "Spdef+ Secretory",
                               "Selenop+/Gstm2high Secretory")))+
  coord_flip()



ggsave(filename = "FIGs2e2_epi_doublet_quant.pdf", plot = epi_doublet, width = 10, height = 16, dpi = 600)





#### Figure Supp 2f: Sample Distribution for All Cells #### 


table <- table(Distal_Named@active.ident ,
               Distal_Named@meta.data$Sample)    # Create a table of counts

df <- data.frame(table) 

table2 <- table(Epi_Named@meta.data$Samlpe)


all_sample_dist <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
                          aes(x = Var1,              # Variable to plot on the x-axis
                              y = Freq,              # Variable to plot on the y-axis
                              fill = factor(Var2,    # Variable to fill the bars
                                            levels = c("mD1","mD2","mD4")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", limits = c("mD1","mD2","mD4"),
                    values = c(natparks.pals(name="Arches2",n=3))) +
  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))               # Text color and horizontal adjustment on y-axis


ggsave(filename = "FIGs2f_all_sample_dist.pdf", plot = epi_sample_dist, width = 16, height = 12, dpi = 600)


#### Figure Supp 2g: Sample Distribution for Epithelial Cells #### 


table <- table(Epi_Named@active.ident ,
               Epi_Named@meta.data$Sample)    # Create a table of counts

df <- data.frame(table) 

table2 <- table(Epi_Named@meta.data$Samlpe)


epi_sample_dist <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
                          aes(x = Var1,              # Variable to plot on the x-axis
                              y = Freq,              # Variable to plot on the y-axis
                              fill = factor(Var2,    # Variable to fill the bars
                                            levels = c("mD1","mD2","mD4")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", limits = c("mD1","mD2","mD4"),
                    values = c(natparks.pals(name="Arches2",n=3))) +
  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))               # Text color and horizontal adjustment on y-axis


ggsave(filename = "FIGs2g_epi_sample_dist.pdf", plot = epi_sample_dist, width = 16, height = 12, dpi = 600)


