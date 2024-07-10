
#### Figure 1: Census of cell types of the mouse uterine tube ####

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

#### Distal and Proximal Datasets ####

Distal <- readRDS(file = "../Distal/20220810_Filtered_Cells.rds" , refhook =  NULL)

Proximal <- readRDS( file = "../Proximal/20220907_Filtered_Cells_Final.rds" , refhook =  NULL)


#### Figure 1b: Cells of the Distal Uterine Tube ####

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

Distal_Named <- SetIdent(Distal_Named, value = Distal_Named@active.ident)



Fibroblasts <- c('#FF9D00' , '#FFB653' , '#FFCB9A')   # Oranges
Muscle <- c('#E55451' , '#FFB7B2') # Reds
Endothelial <- c('#A0E6FF')  # Reds
FiboEpi <- "#FFE0B3" # Reddish Brown
Epi <-c('#6E3E6E','#8A2BE2','#CCCCFF','#DA70D6','#DF73FF','#604791') # Blues/Purples
Immune <- c( '#5A5E6B'  , '#B8C2CC' , '#FC86AA') # Yellowish Brown
Meso <- "#9EFFFF" # Pink
Lut <- "#9DCC00" # Green

colors <- c(Fibroblasts, Muscle, Endothelial, FiboEpi, Epi, Immune, Meso, Lut)


pie(rep(1,length(colors)), col=colors) 


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



p1 <- DimPlot(
  Distal_Named,
  reduction='umap',
  cols=colors,
  pt.size = 0.5,
  label.size = 4,
  label.color = "black",
  repel = TRUE,
  label=F) +
  NoLegend() +
  labs(x="UMAP_1",y="UMAP_2")

LabelClusters(p1, id="ident", color = "black", repel = T , size = 4, box.padding = .75)

ggsave(filename = "FIG1b_all_distal_umap.pdf", plot = p1, width = 8, height = 12, dpi = 600)



## Figure 1c: Distal Uterine Tube Features for Cell Type Identification ##

features <- c("Dcn","Col1a1",           # Fibroblasts        
              "Acta2","Myh11","Myl9",   # Smooth Muscle
              "Pdgfrb","Mcam","Cspg4",  # Pericyte
              "Sele","Vwf","Tek",             # Blood Endothelial
              "Lyve1","Prox1","Icam1",          # Lymph Endothelial
              "Epcam","Krt8",           # Epithelial
              "Foxj1",                  # Ciliated Epithelial
              "Ovgp1",                  # Secretory Epithelial
              "Slc1a3","Pax8","Cd44","Itga6",         # Stem-like Epithelieal 
              "Ptprc",                  # Immune                            
              "Cd8a","Cd4","Cd3e",      # T-Cell                            
              "Klrc1","Runx3",          # T/NK Cell
              "Klrd1",                  # NK Cell
              "Aif1","Cd68","Csf1r","Itgax", # Macrophage
              "Hbb-bs", "Hba-a1",       # Erythrocytes
              "Fras1","Rspo1","Lrrn4","Msln", # Mesothelial
              "Cyp11a1","Bcat1","Fkbp5","Spp1","Prlr") # Luteal Cells

all_dp <- DotPlot(object = Distal_Named,                    # Seurat object
                  assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                  features = features,                  # List of features (select one from above or create a new one)
                  # Colors to be used in the gradient
                  col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
                  col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
                  dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
                  dot.scale = 9,                        # Scale the size of the points
                  group.by = NULL,              # How the cells are going to be grouped
                  split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
                  scale = TRUE,                         # Whether the data is scaled
                  scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
                  scale.min = NA,                       # Set lower limit for scaling
                  scale.max = NA                        # Set upper limit for scaling
)+    
  labs(x = NULL, y = NULL)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  #theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 14 , face = "italic"))+
  theme(axis.text.y = element_text(size = 14))+
  scale_y_discrete(limits = rev(levels(Distal_Named)))+
  theme(legend.title = element_text(size = 14))

ggsave(filename = "FIG1c_all_distal_dotplot.pdf", plot = all_dp, width = 18, height = 10, dpi = 600)

