getwd()
setwd("/Users/Coulter/Desktop/PhD Research/Seurat/Mouse Samples/Distal")
getwd()


#### Load Libraries and Functions ####
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

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(viridisLite)
library(Polychrome)
library(circlize)
library(NatParksPalettes)

#Functions
npcs <- function(seu, var.toal=0.95, reduction="pca"){
  if(is.null(seu@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }
  tmp.var <- (seu@reductions[[reduction]]@stdev)^2
  var.cut <- var.toal*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }
  return(n.pcs)
}

#### Load and Visualize HQ Cells ####

Distal <- readRDS( file = "20220810_Filtered_Cells.rds" , refhook =  NULL)

## Identify Cell Types ##

DoHeatmap(subset(Distal, downsample = 100), 
          features = VariableFeatures(Distal)[1:200], size = 5)+
  scale_fill_gradientn(colors = c( "#a7d5ed", "#e2e2e2", "#e14b31"))+
  theme(text = element_text(size = 4))

DimPlot(object = Distal, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')


DimPlot(object = Distal, 
        reduction = 'umap', 
        group.by = "seurat_clusters",
        cols = c("#EF9245", #1
                 "#EF9245", #2
                 "#FF0000", #3
                 "#EF9245", #4
                 "#FF0000", #5
                 "#EF9245", #6
                 "#FF0000", #7
                 "#FF0000", #8
                 "#FF0000", #9
                 "#EADF75", #10
                 "#DDC700", #11
                 "#FF0000", #12
                 "#EFEF11", #13
                 "#EF9245", #14
                 "#EFEF11", #15
                 "#EFEF11", #16
                 "#EF9245", #17
                 "#EF9245", #18
                 "#FF0000"), #19
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')


features <- c("Dcn","Col1a1",           # Fibroblasts        
              "Acta2","Myh11","Myl9",   # Smooth Muscle
              "Pdgfrb","Mcam","Cspg4",  # Pericyte
              "Sele","Vwf","Tek",             # Blood Endothelial
              "Lyve1","Prox1","Icam1",          # Lymph Endothelial
              "Epcam","Krt8",           # Epithelial
              "Foxj1",                  # Ciliated Epithelial
              "Ovgp1",                  # Secretory Epithelial
              "Slc1a3","Pax8","Krt5","Cd44","Itga6",         # Stem-like Epithelieal 
              "Ptprc",                  # Immune                            
              "Cd8a","Cd4","Cd3e",      # T-Cell                            
              "Klrc1","Runx3",          # T/NK Cell
              "Klrd1","Tnfrsf8","Ncam1",# NK Cell
              "Cd38","Jchain",          # B Cell
              "Aif1","Cd68","Csf1r","Itgax", # Macrophage
              "Hbb-bs", "Hba-a1",       # Erythrocytes
              "Cyp11a1","Bcat1","Fkbp5","Spp1","Prlr", # Luteal Cells
              "Fras1","Rspo1","Lrrn4","Msln", # Mesothelial
              "Crisp3","Bche")

DotPlot(object = Distal,                    # Seurat object
        assay = 'RNA',                        # Name of assay to use.  Default is the active assay
        features = "Krt5",                  # List of features (select one from above or create a new one)
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
  theme(legend.title = element_text(size = 12))

## Rename ##

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
                               '18' = "Lymph Endothelial/Epithelial")

## Refined Feature Plot ##

features <- c("Dcn","Col1a1",           # Fibroblasts        
              "Acta2","Myh11","Myl9",   # Smooth Muscle
              "Pdgfrb","Mcam","Cspg4",  # Pericyte
              "Sele","Vwf",             # Blood Endothelial
              "Lyve1","Prox1",          # Lymph Endothelial
              "Epcam","Krt8",           # Epithelial
              "Foxj1",                  # Ciliated Epithelial
              "Ovgp1",                  # Secretory Epithelial
              "Slc1a3","Pax8","Krt5", "Cd44","Itga6",          # Stem-like Epithelieal 
              "Fras1","Rspo1","Lrrn4","Msln", # Mesothelial
              "Ptprc",                  # Immune                            
              "Cd8a","Cd4","Cd3e",      # T-Cell                            
              "Klrc1","Runx3",          # T/NK Cell
              "Aif1","Cd68","Csf1r","Itgax", # Macrophage
              "Hbb-bs", "Hba-a1",       # Erythrocytes
              "Cyp11a1","Bcat1","Fkbp5") # Luteal Cells



dp <- DotPlot(object = Distal_Named,                    # Seurat object
        assay = 'RNA',                        # Name of assay to use.  Default is the active assay
        features = features,                  # List of features (select one from above or create a new one)
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
           y = NULL) +                            # y-axis label
  RotatedAxis()+
  scale_y_discrete(limits = c("Luteal",
                              "Erythrocyte",
                              "Macrophage",
                              "T/NK Cell",
                              "Mesothelial",
                              "Intermediate Epithelial",
                              "Stem-like Epithelial 2",
                              "Stem-like Epithelial 1",
                              "Secretory Epithelial",
                              "Ciliated Epithelial 2",
                              "Ciliated Epithelial 1",
                              "Lymph Endothelial/Epithelial",
                              "Blood Endothelial",
                              "Pericyte",
                              "Smooth Muscle",
                              "Epithelial/Fibroblast",
                              "Fibroblast 3",
                              "Fibroblast 2",
                              "Fibroblast 1"))+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 12 , face = "italic"))+
  theme(axis.text.y = element_text(size = 12))+
  theme(legend.title = element_text(size = 12))

ggsave(filename = "all_distal_dot_plot.pdf", plot = dp, width = 15, height = 8, dpi = 600)

## Stacked Bar Graph ##

table <- table(Distal_Named@active.ident ,
               Distal_Named@meta.data$Sample)    # Create a table of counts

df <- data.frame(table) 

table2 <- table(Distal_Named@meta.data$Samlpe)


sample_dist <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
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
  # Change the order of the x-axis
  # scale_x_discrete(limits = c("FAPs", "Tenocytes", "Endothelial cells", "Smooth muscle cells", "Mature skeletal muscle", "MuSCs and progenitors", 
  #  "Schwann and Neural/Glial cells", "Pro-inflammatory macrophages", "Anti-inflammatory macrophages", "Monocyte (Cxcl10+)", 
  #  "Monocyte (Patrolling)", "Macrophage (Self-renewing)", "Neutrophils", "Dendritic cells", "B cells", "T cells", "NK cells")) +
  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1)) +              # Text color and horizontal adjustment on y-axis
  # Split the plot by Age
  scale_x_discrete(limits = c( "Fibroblast 1", 
                               "Fibroblast 2", 
                               "Fibroblast 3",
                               "Epithelial/Fibroblast",
                               "Smooth Muscle",
                               "Pericyte",
                               "Blood Endothelial",
                               "Ciliated Epithelial 1",
                               "Ciliated Epithelial 2",
                               "Secretory Epithelial",
                               "Stem-like Epithelial 1",
                               "Stem-like Epithelial 2",
                               "Intermediate Epithelial",
                               "Mesothelial",
                               "T/NK Cell",
                               "Macrophage",
                               "Erythrocyte",
                               "Luteal"))
ggsave(filename = "All_sample_dist.pdf", plot = sample_dist, width = 16, height = 12, dpi = 600)

# Quantify #
write.csv(df, file = "quantify_distal.csv")


## UMAP Plot ##

DimPlot(object = Distal_Named,                # Seurat object  
        reduction = 'umap',                 # Axes for the plot (UMAP, PCA, etc.) 
        #group.by = "Patient",       # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
        repel = TRUE,                       # Whether to repel the cluster labels
        label = FALSE,                       # Whether to have cluster labels 
        cols = c("#47816B", #1
                 "#889771", #2
                 "#46BDC6", #3
                 "#D5624E", #4
                 "#6938E0", #5
                 "#C7C577", #6
                 "#6D8CE2", #7
                 "#3364DD", #8
                 "#452A9B", #9
                 "#CE6D5C", #10
                 "#943E49", #11
                 "#A289DB", #12
                 "#E79624", #13
                 "#9BB9A1", #14
                 "#D6BA71", #15
                 "#EDA194", #16
                 "#E256C5", #17
                 "#66E6ED", #18
                 "#010A0A"), #19
        pt.size = 0.4,                      # Size of each dot is (0.1 is the smallest)
        label.size = 0.5) +                   # Font size for labels    
  # You can add any ggplot2 1customizations here
  labs(title = 'Colored by Cluster')        # Plot title
#### Subset Epithelial Cells ####

FeaturePlot(Distal_Named, features = c("Epcam"), pt.size = 0.75)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

Epi_Filter <- subset(Distal_Named, 
                     idents = c("Ciliated Epithelial 1",
                                "Ciliated Epithelial 2",
                                "Secretory Epithelial",
                                "Stem-like Epithelial 1",
                                "Stem-like Epithelial 2",
                                "Intermediate Epithelial"))

# Need to Normalize and Run PCA before Running Harmony ##

Epi_Filter <- NormalizeData(object = Epi_Filter, assay = 'RNA')
Epi_Filter <- FindVariableFeatures(object = Epi_Filter, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
Epi_Filter <- ScaleData(object = Epi_Filter, assay = 'RNA')
Epi_Filter <- RunPCA(object = Epi_Filter, assay = 'RNA', features = VariableFeatures(object = Epi_Filter),
                     reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')

## Determine the 'dimensionality' of the dataset ##

ElbowPlot(Epi_Filter,
          reduction = 'pca_RNA',
          ndims = 50) 

# Calculate the number of PCs that contain some proportion (95%) of the variance
npcs <- function(seu, var.toal=0.95, reduction="pca"){
  if(is.null(seu@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }
  tmp.var <- (seu@reductions[[reduction]]@stdev)^2
  var.cut <- var.toal*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }
  return(n.pcs)
}

pcs <- npcs(seu = Epi_Filter, var.toal = 0.95, reduction = 'pca_RNA')

## Run Harmony ##

gc()

Epi_Filter <- RunHarmony(Epi_Filter, 
                         group.by.vars = c("Sample"), 
                         reduction = "pca_RNA", 
                         assay = "RNA", 
                         plot_convergence = TRUE)

## Normal Seurat clustering ##

Epi_Filter <- FindNeighbors(object = Epi_Filter, 
                            reduction = "harmony",
                            k.param = 50, 
                            dims = 1:pcs, 
                            graph.name = 'harmony_snn')

Epi_Filter <- FindClusters(object = Epi_Filter, 
                           resolution = 0.7, 
                           graph.name = 'harmony_snn')


## Set Python UMAP via reticulate ##

umap.method = 'umap-learn'
metric = 'correlation'

Epi_Filter <- RunUMAP(object = Epi_Filter, 
                      reduction = "harmony", dims = 1:pcs)

## UMAP ##

DimPlot(object = Epi_Filter, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

FeaturePlot(Epi_Filter, features = c("Krt5"), pt.size = 0.4)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

## Save as RDS file ##

saveRDS(Epi_Filter, file = "20220817_Distal_Epi_Cells.rds")

#### Epithelial Analysis ####

Epi_Filter <- readRDS(file = "20220817_Distal_Epi_Cells.rds" , refhook =  NULL)

EpiAllMarkers <- FindAllMarkers(Epi_Filter)
write.csv(EpiAllMarkers, file = "20220817_Distal_Epi_All_Markers.csv")

features <- c("Krt8","Epcam",
              "Foxj1","Fam183b",
              "Dnali1",
              "Ovgp1","Pax8","Sox17",
              "Egr1",
              "Msln","Acta2","Thy1","Timp3","Lgals1","Dcn",
              "Slc1a3","Cd44","Sox9","Krt5","Krt14","Krt17", "Mki67")

stained_features <- c("Cd44","Cdkn2a","Ctnnb1","Itga6","Krt5","Krt6b","Krt6a",
                      "Krt7","Lef1","Mki67","Pax2","Pax8","Prom1","Tp53","Wt1")

stem_features <- c("Krt5","Krt17","Cd44","Prom1","Kit","Aldh1a1","Aldh1a2","Aldh1a3",
                   "Efnb1","Ephb1","Trp63","Sox2","Sox9","Klf4","Rnf43","Foxm1",
                   "Pax8","Nanog","Itga6","Psca","Tcf3","Tcf4","Nrp1","Slc1a3","Tnfrsf19",
                   "Smo","Lrig1","Ezh2","Egr1","Tacstd2","Dusp1","Slc38a2","Malat1",
                   "Btg2","Cdkn1c","Pdk4","Nedd9","Fos","Jun","Junb","Zfp36",
                   "Neat1","Gadd45g","Gadd45b")

## Visualize Distribution of Expressed Genes ##

DotPlot(object = Epi_Filter,                    # Seurat object
        assay = 'RNA',                        # Name of assay to use.  Default is the active assay
        features = features,                  # List of features (select one from above or create a new one)
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
  theme(legend.title = element_text(size = 12))


VlnPlot(Epi_Filter, features = "Stat3")

## Rename ##

Epi_Named <- RenameIdents(Epi_Filter, 
                             '0' = "Spdef+ Secretory", 
                             '1' = "Slc1a3+ Progenitor", 
                             '2' = "Msln+ Progenitor",
                             '3' = "Ciliated 1", 
                             '4' = "Ciliated 2", 
                             '5' = "Pax8low/Prom1+ Cilia-forming", 
                             '6' = "Fibroblast-like",
                             '7' = "Slc1a3med/Sox9+ Cilia-forming",
                             '8' = "Selenop+/Gstm2high Secretory")

Epi_Named@active.ident <- factor(x = Epi_Named@active.ident, levels = c( c("Slc1a3+ Progenitor",
                                                             "Msln+ Progenitor",
                                                             "Slc1a3med/Sox9+ Cilia-forming",
                                                             "Pax8low/Prom1+ Cilia-forming", 
                                                             "Fibroblast-like",
                                                             "Spdef+ Secretory",
                                                             "Selenop+/Gstm2high Secretory",
                                                             "Ciliated 1",
                                                             "Ciliated 2")))

VlnPlot(Epi_Named, 
        features = "Krt17",
        slot = "data",
        pt.size = 1,
        cols = c("#35EFEF", #1
                 "#00A1C6", #2
                 "#2188F7", #3
                 "#EA68E1", #4
                 "#59D1AF", #5
                 "#B20224", #6
                 "#F28D86", #7
                 "#A374B5", #8
                 "#9000C6"))+ #9
  theme(plot.title = element_text(size = 32, face = "bold.italic"))+
  scale_x_discrete(limits = c("Slc1a3+ Progenitor",
                              "Msln+ Progenitor",
                              "Slc1a3med/Sox9+ Cilia-forming",
                              "Pax8low/Prom1+ Cilia-forming", 
                              "Fibroblast-like",
                              "Spdef+ Secretory",
                              "Selenop+/Gstm2high Secretory",
                              "Ciliated 1",
                              "Ciliated 2"))+
  theme(axis.text.x = element_text(size = 16, angle = 60))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.y.left = element_text(size = 16))
  
ggsave(filename = "Foxj1_epi_umap.pdf", plot = fj, width = 15, height = 10, dpi = 600)

## Analyze Named Clusters ##

named_features <- c("Krt8","Epcam",
                    "Slc1a3","Cd44","Sox9",
                    "Ovgp1","Sox17","Pax8", "Egr1",
                    "Itga6", "Bmpr1b",
                    "Msln", "Anpep", "Klf6",
                    "Dpp6", "Sec14l3", "Fam161a",
                    "Prom1", "Ly6a", "Kctd8", "Adam8",
                    "Foxj1","Rgs22","Fam183b",
                    "Dnali1", "Mt1" , "Dynlrb2",
                    "Dcn", "Col1a1", "Col1a2", "Timp3", "Pdgfra","Lgals1",
                    "Upk1a", "Thrsp",
                    "Selenop", "Gstm2",
                    "Krt5","Krt14","Krt17")

Wnt_targets <- c('Lgr4' , 'Enc1' , 'Bcl2l2' , 'Yap1',
                 'Lef1' , 'Fzd7' , 'Axin2' , 'Ccnd1' , 'Dkk1' , 'Lgr5',
                 'Tgif1' , 'Myc' , 'Msln' , 'Postn' , 'Procr' , 'Mmp3' , 'Wnt1',
                 'Rhou','Kif3a')


miR34 <- c('Sirt1' , 'Bcl2' , 'Ccnd1' , 'Ccne2' , 'Cdk4' , 'Cdk6' , 'Myc' , 
           'E2f3' , 'Met' , 'Mycn',
           'Rhoa', 'Ccp110', 'Rras',
           'Arhgap1' , 'Arhgdib','Mir34a','Mir34b','Mir34c')


DotPlot(object = Epi_Named,                    # Seurat object
        assay = 'RNA',                        # Name of assay to use.  Default is the active assay
        features = "Krt5",                  # List of features (select one from above or create a new one)
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
  scale_y_discrete(limits = rev(c("Slc1a3+ Progenitor",
                     "Msln+ Progenitor",
                     "Slc1a3med/Sox9+ Cilia-forming",
                     "Pax8low/Prom1+ Cilia-forming",
                     "Ciliated 1",
                     "Ciliated 2", 
                     "Fibroblast-like",
                     "Spdef+ Secretory",
                     "Selenop+/Gstm2high Secretory")))+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 16 , face = "italic"))+
  theme(axis.text.y = element_text(size = 16))+
  theme(legend.title = element_text(size = 16))+ RotatedAxis()


RidgePlot(Epi_Named, features = "Krt5", ncol = 2)

DimPlot(object = Epi_Named,                # Seurat object  
        reduction = 'umap',                 # Axes for the plot (UMAP, PCA, etc.) 
        #group.by = "Patient",       # Labels to color the cells by ("seurat_clusters", "Age", "Time.Point)  
        repel = TRUE,                       # Whether to repel the cluster labels
        label = FALSE,                       # Whether to have cluster labels 
        cols = c("#B20224", #1
                 "#35EFEF", #2
                 "#00A1C6", #3
                 "#A374B5", #4
                 "#9000C6", #5
                 "#EA68E1", #6
                 "#59D1AF", #7
                 "#2188F7", #8
                 "#F28D86"), #9
                 
        pt.size = 0.6,                      # Size of each dot is (0.1 is the smallest)
        label.size = 0.5) +                   # Font size for labels    
  # You can add any ggplot2 1customizations here
  labs(title = 'Colored by Cluster')        # Plot title


## Stacked Bar Graph ##

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
  # Change the order of the x-axis
  # scale_x_discrete(limits = c("FAPs", "Tenocytes", "Endothelial cells", "Smooth muscle cells", "Mature skeletal muscle", "MuSCs and progenitors", 
  #  "Schwann and Neural/Glial cells", "Pro-inflammatory macrophages", "Anti-inflammatory macrophages", "Monocyte (Cxcl10+)", 
  #  "Monocyte (Patrolling)", "Macrophage (Self-renewing)", "Neutrophils", "Dendritic cells", "B cells", "T cells", "NK cells")) +
  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))               # Text color and horizontal adjustment on y-axis


ggsave(filename = "All_epi_sample_dist.pdf", plot = epi_sample_dist, width = 16, height = 12, dpi = 600)
  
write.csv( df, file = "distal_epi_distribution.csv")



#### PHATE Analysis ####

#https://github.com/KrishnaswamyLab/PHATE/issues/104

Distal_PHATE <- RunPHATE(Epi_Named, knn = 80, reduction = "harmony", npca = 39)

DimPlot(Distal_PHATE , reduction = "phate", 
        cols = c("#B20224", #1
                 "#35EFEF", #2
                 "#00A1C6", #3
                 "#A374B5", #4
                 "#9000C6", #5
                 "#EA68E1", #6
                 "#59D1AF", #7
                 "#2188F7", #8
                 "#F28D86"),
        pt.size = 0.7,
        shuffle = TRUE,
        seed = 0,
        label = FALSE)

saveRDS(Epi_Filter, file = "20220826_Distal_Epi_PHATE.rds")

FeaturePlot(Distal_PHATE, features = "Epcam", reduction = "phate")

Beeg_PHATE <- Distal_PHATE

Beeg_PHATE@reductions[["phate"]]@cell.embeddings <- Distal_PHATE@reductions[["phate"]]@cell.embeddings*100

FeaturePlot(Beeg_PHATE, features = c("Krt17"), reduction = "phate", pt.size = 0.9, order = T)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

#### Monocle3 Analysis ####
# Prepare cds File #


library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(tidyverse)

library(BiocManager)
library(remotes)
library(devtools)


gene_annotation <- as.data.frame(rownames(Distal_PHATE@reductions[["harmony"]]@feature.loadings),
                                 row.names = rownames(Distal_PHATE@reductions[["harmony"]]@feature.loadings))

colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(Distal_PHATE@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = Distal_PHATE@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "sample"


cell_metadata <- separate(cell_metadata, col = 'sample', 
                          into = c('Patient', 'Barcode'), sep = '_')

New_matrix <- Distal_PHATE@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(Distal_PHATE@reductions[["harmony"]]@feature.loadings), ]



expression_matrix <- New_matrix


cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


list_cluster <- Distal_PHATE@active.ident
names(list_cluster) <- Distal_PHATE@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <- (Distal_PHATE@reductions[["phate"]]@cell.embeddings*100)

cds_from_seurat

# Process cds file #

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4,
           group_label_size = 0,
           cell_size = 0.9,)+ 
  scale_color_manual(values = c("#B20224", #1
                                "#35EFEF", #2
                                "#00A1C6", #3
                                "#A374B5", #4
                                "#9000C6", #5
                                "#EA68E1", #6
                                "#59D1AF", #7
                                "#2188F7", #8
                                "#F28D86"))

saveRDS(cds_from_seurat, file = "20221101_Distal_Epi_PHATE_Monocle3.rds")

## Visualize ##

cds <- order_cells(cds_from_seurat)

saveRDS(cds, file = "20221101_Distal_Epi_PHATE_Monocle3.rds")

cds <- readRDS(file = "20221101_Distal_Epi_PHATE_Monocle3.rds" , refhook = NULL)

x <- pseudotime(cds)

y <- as.data.frame(x)


gene_fits <- fit_models(cds, model_formula_str = "~embryo.time")


z <- cds@clusters$UMAP$clusters

zz <- as.data.frame(z)

total <- rbind(y, zz)


a = merge(y, zz)


plot_genes_violin(cds, group_cells_by=pseudotime(cds), ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


plot_cells(cds_from_seurat,
           genes="Dcn",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           cell_size = .01,
           cell_stroke = 1)+
  scale_color_viridis_c(option="F",begin=.2,end=0.9, direction = -1)

cds <- order_cells(cds_from_seurat)

plot_cells(cds, 
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           cell_size = .01,
           cell_stroke = 1)



plot_cells(cds,
           genes="Slc1a3",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           cell_size = .01,
           cell_stroke = 2)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)


cds_sub <- choose_graph_segments(cds)

pr_graph_test_res <- graph_test(cds_sub, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))



cds_sub <- preprocess_cds(cds_sub, num_dim = 39)
cds_sub <- reduce_dimension(cds_sub)
cds_sub <- cluster_cells(cds_sub)
cds_sub <- learn_graph(cds_sub, use_partition = F)


plot_cells(cds_sub, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4,
           group_label_size = 0,
           cell_size = 0.9,)

plot_cells(cds_sub,
           genes="Krt5",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           cell_size = .01,
           cell_stroke = 2)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)

















#### Cell Cycle Analysis ####

library(biomaRt)
library(httr)

## Convert Mouse to Human ##
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)



epi_cycle <- CellCycleScoring(Epi_Named, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)

## Save as RDS file ##

saveRDS(epi_cycle, file = "20220818_Distal_Epi_Cell_Cycle.rds")

epi_cycle <- readRDS(file = "20220818_Distal_Epi_Cell_Cycle.rds" , refhook = NULL)

### Analysis ###

VlnPlot(epi_cycle, features = c("Ect2", "Cdk4", "Mcm6","Mki67"), 
           split.by = "Phase",
           group.by = "seurat_clusters", 
           ncol = 2, slot = "counts")

VlnPlot(epi_cycle, features = c("Pax8"), 
        split.by = "Phase",
        group.by = "seurat_clusters", slot = "counts")

table(epi_cycle$Phase)

DimPlot(epi_cycle)


## Stacked Bar ##


table <- table(epi_cycle@meta.data$old.ident ,
               epi_cycle@meta.data$Phase)    # Create a table of counts

df <- data.frame(table) 

ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
       aes(x = Var1,              # Variable to plot on the x-axis
           y = Freq,              # Variable to plot on the y-axis
           fill = factor(Var2,    # Variable to fill the bars
                         levels = c("G1","G2M","S")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", limits = c("G1","G2M","S"),
                    values = c(natparks.pals(name="Arches2",n=3))) +
  # Change the order of the x-axis
  # scale_x_discrete(limits = c("FAPs", "Tenocytes", "Endothelial cells", "Smooth muscle cells", "Mature skeletal muscle", "MuSCs and progenitors", 
  #  "Schwann and Neural/Glial cells", "Pro-inflammatory macrophages", "Anti-inflammatory macrophages", "Monocyte (Cxcl10+)", 
  #  "Monocyte (Patrolling)", "Macrophage (Self-renewing)", "Neutrophils", "Dendritic cells", "B cells", "T cells", "NK cells")) +
  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))+              # Text color and horizontal adjustment on y-axis
  scale_x_discrete(limits = c ("Progenitor",
                               "Msln+ Progenitor",
                               "Slc1a3+/Sox9+ Cilia-forming",
                               "Pax8+/Prom1+ Cilia-forming", 
                               "Fibroblast-like",
                               "Secretory 1",
                               "Secretory 2",
                               "Ciliated 1",
                               "Ciliated 2"))             # Text color and horizontal adjustment on y-axis



# Split the plot by Age
# scale_x_discrete(limits = c( "Fibroblast 1", 
"Fibroblast 2", 
"Fibroblast 3",
"Epithelial/Fibroblast",
"Smooth Muscle",
"Pericyte",
"Blood Endothelial",
"Lymph Endothelial/Epithelial",
"Ciliated Epithelial 1",
"Ciliated Epithelial 2",
"Secretory Epithelial",
"Stem-like Epithelial 1",
"Stem-like Epithelial 2",
"Intermediate Epithelial",
"Mesothelial",
"T/NK Cell",
"Macrophage",
"Erythrocyte",
"Luteal"))



# view cell cycle scores and phase assignments
head(epi_cycle[[]])


I just use an R function posted here:
  https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
  
  and then go
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)


## Doublet Analysis ##

# Doublet Dataframe #

pANN <- coalesce(Epi_Filter@meta.data$pANN_0.25_0.005_191, 
                 Epi_Filter@meta.data$pANN_0.25_0.3_101,
                 Epi_Filter@meta.data$pANN_0.25_0.3_133)

doublet <- Epi_Filter@meta.data$Doublet

cluster <-  Epi_Filter@meta.data$harmony_snn_res.0.7

pANN_df <- data.frame(cluster,pANN)

Avg_pANN <- pANN_df %>%
  group_by(cluster) %>%
  summarise_at(vars(pANN), list(name = mean))

doublet_table <- table(cluster,doublet)
doublet_df <- data.frame(doublet_table)

Pct_doublet <- transform(doublet_df,            # Calculate percentage by group
                         perc = ave(Freq,
                                    cluster,
                                    FUN = prop.table))

Avg_pANN

write.csv(Avg_pANN , "Distal_Avg_pANN.csv")
write.csv(Pct_doublet , "Distal_Pct_doublet.csv")

