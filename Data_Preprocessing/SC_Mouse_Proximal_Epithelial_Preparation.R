getwd()
setwd("/Users/Coulter/Desktop/PhD Research/Seurat/Mouse Samples/Proximal")
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
library(DropletUtils)

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

Proximal <- readRDS( file = "20220907_Filtered_Cells_Final.rds" , refhook =  NULL)

## Identify Cell Types ##

DoHeatmap(subset(Proximal, downsample = 100), 
          features = VariableFeatures(Proximal)[1:200], size = 5)+
  scale_fill_gradientn(colors = c( "#a7d5ed", "#e2e2e2", "#e14b31"))+
  theme(text = element_text(size = 4))

DimPlot(object = Proximal, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = Proximal, 
        reduction = 'umap', 
        group.by = "seurat_clusters",
        cols = c(  "#09CE45", #1
                   "#005CFF", #2
                   "#09CE45", #3
                   "#09CE45", #4
                   "#9C3CAA", #5
                   "#303AF2", #6
                   "#09CE45", #7
                   "#005CFF", #8
                   "#09CE45", #9
                   "#09CE45"), #10
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.9, 
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
              "Fras1","Rspo1","Lrrn4","Msln") # Mesothelial

DotPlot(object = Proximal,                    # Seurat object
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

## Rename ##

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



DotPlot(object = Proximal_Named,                    # Seurat object
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
  scale_y_discrete(limits = c("Immune",
                              "Mesothelial",
                              "Stem-like Epithelial",
                              "Secretory Epithelial",
                              "Ciliated Epithelial",
                              "Endothelial",
                              "Smooth Muscle",
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

## Stacked Bar Graph ##

table <- table(Proximal_Named@active.ident ,
               Proximal_Named@meta.data$Sample)    # Create a table of counts

df <- data.frame(table) 

table2 <- table(Proximal_Named@meta.data$Samlpe)


ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
       aes(x = Var1,              # Variable to plot on the x-axis
           y = Freq,              # Variable to plot on the y-axis
           fill = factor(Var2,    # Variable to fill the bars
                         levels = c("mP2","mP3","mP4")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", limits = c("mP2","mP3","mP4"),
                    values = c(natparks.pals(name="Glacier",n=3))) +
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
                               "Smooth Muscle",
                               "Endothelial",
                               "Ciliated Epithelial",
                               "Secretory Epithelial",
                               "Stem-like Epithelial",
                               "Mesothelial",
                               "Immune"))

# Quantify #
write.csv(df, file = "quantify_proximal.csv")


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

FeaturePlot(Proximal_Named, features = c("Epcam"), pt.size = 0.75)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

Epi_Filter <- subset(Proximal_Named, 
                     idents = c("Ciliated Epithelial",
                                "Secretory Epithelial",
                                "Stem-like Epithelial"))

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
                            k.param = 30, 
                            dims = 1:pcs, 
                            graph.name = 'harmony_snn')

Epi_Filter <- FindClusters(object = Epi_Filter, 
                           resolution = 0.9, 
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
        pt.size = 1.3, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

FeaturePlot(Epi_Filter, features = c("Slc1a3"), pt.size = 0.4)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

## Save as RDS file ##

saveRDS(Epi_Filter, file = "20220914_Proximal_Epi_Cells.rds")

#### Epithelial Analysis ####

Epi_Filter <- readRDS(file = "20220914_Proximal_Epi_Cells.rds" , refhook =  NULL)

EpiAllMarkers <- FindAllMarkers(Epi_Filter)
write.csv(EpiAllMarkers, file = "20220914_Proximal_Epi_All_Markers_2.csv")

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

distal_features <- c("Krt8","Epcam",
                    "Slc1a3","Cd44","Sox9",
                    "Ovgp1","Sox17","Pax8", "Egr1",
                    "Itga6", "Bmpr1b",
                    "Msln", "Anpep", "Klf6",
                    "Dpp6", "Sec14l3", "Fam161a",
                    "Prom1", "Ly6a", "Kctd8", "Adam8",
                    "Dcn", "Col1a1", "Col1a2", "Timp3", "Pdgfra","Lgals1",
                    "Upk1a", "Thrsp",
                    "Selenop", "Gstm2",
                    "Foxj1","Fam183b",
                    "Rgs22","Dnali1", "Mt1" , "Dynlrb2",
                    "Krt5","Krt14","Krt17")

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


VlnPlot(Epi_Filter, features = "Gsto1")

## Rename ##

Epi_Named <- RenameIdents(Epi_Filter, 
                          '0' = "Spdef+ Secretory", 
                          '1' = "Slc1a3+ Progenitor", 
                          '2' = "Itga6+ Secretory",
                          '3' = "Pax8+/Prom1+ Cilia-forming", 
                          '4' = "Sox17high Secretory", 
                          '5' = "Kcne3+ Secretory")

Epi_Named@active.ident <- factor(x = Epi_Named@active.ident, levels = c("Slc1a3+ Progenitor",
                                                                           "Pax8+/Prom1+ Cilia-forming",
                                                                           "Spdef+ Secretory",
                                                                           "Itga6+ Secretory",
                                                                           "Sox17high Secretory",
                                                                           "Kcne3+ Secretory"))

VlnPlot(Epi_Named, 
        features = "Foxa2",
        slot = "data",
        pt.size = 0,
        cols = c("#35EFEF", #1
                 "#E95FE0", #2
                 "#B20224", #3
                 "#F28D86", #4
                 "#FB1111", #5
                 "#FEB0DB"))+ #6
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  scale_x_discrete(limits = c ("Slc1a3+ Progenitor",
                               "Pax8+/Prom1+ Cilia-forming",
                               "Dbi+/Spdefhigh Secretory",
                               "Wfdc2+ Secretory",
                               "Sox17high Secretory",
                               "Kcne3+ Secretory"))



## Analyze Named Clusters ##

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
                    "Rgs22","Dnali1", "Mt1" , "Dynlrb2",
                    "Krt5","Krt14","Krt17")


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
  scale_y_discrete(limits = c("Pax8+/Prom1+ Cilia-forming",
                              "Kcne3+ Secretory",
                               "Sox17high Secretory",
                               "Spdef+ Secretory",
                               "Itga6+ Secretory",
                               "Slc1a3+ Progenitor"
                               ))

ggsave(filename = "20231027_epi_proximal_dotplot.pdf", plot = prox_dp, width = 18, height = 10, dpi = 600)


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
  labs(title = 'Colored by Cluster')        # Plot title

ggsave(filename = "epi_proximal_umap.pdf", plot = epi_umap, width = 15, height = 12, dpi = 600)



## Stacked Bar Graph ##

table <- table(Epi_Named@active.ident ,
               Epi_Named@meta.data$Sample)    # Create a table of counts

df <- data.frame(table) 

table2 <- table(Epi_Named@meta.data$Samlpe)


ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
       aes(x = Var1,              # Variable to plot on the x-axis
           y = Freq,              # Variable to plot on the y-axis
           fill = factor(Var2,    # Variable to fill the bars
                         levels = c("mP2","mP3","mP4")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", limits = c("mP2","mP3","mP4"),
                    values = c(natparks.pals(name="Glacier",n=3))) +
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
  scale_x_discrete(limits = c ("Progenitor",
                               "Secretory 1",
                               "Secretory 2",
                               "Secretory 3",
                               "Secretory 4",
                               "Pax8+/Prom1+ Cilia-forming"))

write.csv( df, file = "distal_epi_distribution.csv")




#### PHATE ####

Proximal_PHATE <- RunPHATE(Epi_Named, knn = 65, reduction = "harmony", npca = 45)

DimPlot(Proximal_PHATE , reduction = "phate", 
        cols = c("#B20224", #1
                 "#35EFEF", #2
                 "#F28D86", #3
                 "#E95FE0", #4
                 "#FB1111", #5
                 "#FEB0DB"),
        pt.size = 2.4,
        shuffle = TRUE,
        seed = 0,
        label = FALSE)

FeaturePlot(Proximal_PHATE, features = "Epcam", reduction = "phate")

Beeg_PHATE <- Proximal_PHATE

Beeg_PHATE@reductions[["phate"]]@cell.embeddings <- Proximal_PHATE@reductions[["phate"]]@cell.embeddings*100

FeaturePlot(Beeg_PHATE, features = c("Prom1"), reduction = "phate", pt.size = 2.4)+
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


gene_annotation <- as.data.frame(rownames(Proximal_PHATE@reductions[["harmony"]]@feature.loadings),
                                 row.names = rownames(Proximal_PHATE@reductions[["harmony"]]@feature.loadings))

colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(Proximal_PHATE@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = Proximal_PHATE@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "sample"


cell_metadata <- separate(cell_metadata, col = 'sample', 
                          into = c('Patient', 'Barcode'), sep = '_')

New_matrix <- Proximal_PHATE@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(Proximal_PHATE@reductions[["harmony"]]@feature.loadings), ]



expression_matrix <- New_matrix


cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


list_cluster <- Proximal_PHATE@active.ident
names(list_cluster) <- Proximal_PHATE@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <- (Proximal_PHATE@reductions[["phate"]]@cell.embeddings*100)

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
           cell_size = 2.4,)+ 
  scale_color_manual(values = c("#B20224", #1
                                "#35EFEF", #2
                                "#F28D86", #3
                                "#E95FE0", #4
                                "#FB1111", #5
                                "#FEB0DB"))

saveRDS(Proximal_PHATE, file = "20220914_Proximal_Epi_PHATE.rds")

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

plot_cells(cds,
           genes="Foxj1",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           cell_size = .01,
           cell_stroke = 2)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)





