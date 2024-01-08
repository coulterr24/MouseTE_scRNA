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

## Preprocessing and Filtering ##

#### Read in Metadata ####

meta <- fread("Mouse_Proximal_Metadata.csv")

#### SoupX Preprocessing ####

soup.list <- list()          # Prepare Empty List         
soup.mat.list <- list()      # Prepare Empty List

## Load souplist objects from cellranger outputs ##

soup.list <- lapply( 
  as.list(paste0(meta$data.dir)),  # Load in Folders with 'Filtered','Raw', and 'Analysis'
  FUN = SoupX::load10X,            # Loads and  Identifies which droplets are cells
  keepDroplets=TRUE
)                                     

## Estimate the Background Contamination ##

soup.list.est <- lapply(
  soup.list,
  FUN = function(sc){
    return(tryCatch(autoEstCont(sc), error=function(e) NULL))
  }
)

## Remove Background to Calculate the Resulting Corrected Count Matrix ##

adj.mat.list <- lapply(
  soup.list.est,
  FUN = function(sc){
    return(tryCatch(adjustCounts(sc), error=function(e) NULL))
  }
)

## Save Adjust Matrices as New Folder ## 

for(i in 1:length(adj.mat.list)){
  if(!is.null(adj.mat.list[[i]]) & !file.exists(paste0(meta$data.dir[i],"/soupx/matrix.mtx"))){
    DropletUtils:::write10xCounts(
      path=paste0(meta$data.dir[i],"/soupx"), # Path to each sample's cellranger count output
      adj.mat.list[[i]]
    )
  }
}

## save Rho values ##

rhos <- list() # 'rho' is parametrized as the 'contamination fraction' 

for(i in 1:length(soup.list.est)){
  rhos[[i]] <- mean(soup.list.est[[i]]$metaData$rho) # Append Metadata w/ Rhos
}
rhos <- do.call(rbind,rhos)
meta$soupx.rho <- rhos

## Save metadata file with the Rho values ##

write.table(x = meta, 
            file = "Mouse_Proximal_Metadata_soup.csv",
            sep = ",",
            row.names = FALSE)

#### Load in Count Matrices ####

mat.list <- list()      # Prepare Empty List 
soupx.used <- list()    # Prepare Empty List 

## Read 10X Data After SoupX ##

for(i in 1:length(meta$data.dir)){ 
  if(file.exists(paste0(meta$data.dir[i], '/soupx'))){ #first try soupx
    cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
    mat.list[[i]] <- Read10X(data.dir = paste0(meta$data.dir[i],"/soupx"))
    soupx.used[[i]] <- T
  }else if(file.exists(paste0(meta$data.dir[i], '/filtered_feature_bc_matrix'))){ #then try cellranger outputs
    cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
    mat.list[[i]] <- Read10X(data.dir = paste0(meta$data.dir[i], '/filtered_feature_bc_matrix'))
    soupx.used[[i]] <- F
  }else{
    cat("Data not found for # ", i, " (", meta$data.dir[i], ")", "\n")
    soupx.used[[i]] <- NULL
  }
}

meta$soupx.used <- unlist(soupx.used) # Checks if SoupX was used
cat(sum(unlist(lapply(mat.list, ncol))),"cells (total) loaded...\n") # Call the number of cells loaded

#### Basic Seurat Workflow ####

seur.list <- list() # Prepare Empty List 

## Initialize Seurat Objects ##

seur.list <- lapply( 
  mat.list,
  FUN = function(mat){
    return(CreateSeuratObject(
      counts = mat, 
      project = 'sc-mouse'
    ))
  }
)  

## Modify Seurat Objects ##

# Read Table of Mt-genes #

Mt_Table <- read.csv("/Users/Coulter/Desktop/PhD Research/Seurat/Mouse Samples/mt_mouse_genes.csv") # Table of mt-genes

mt_df <- data.frame(unique(Mt_Table$TrnP)) # Unique Genes in a List

mt_gene_list <- mt_df$unique.Mt_Table.TrnP.

for(i in 1:length(seur.list)){
  cat(' #####################################\n',
      '### Processing dataset number ', i, '###\n',
      '#####################################\n')
  # Add meta data
  for(md in colnames(meta)){
    seur.list[[i]][[md]] <- meta[[md]][i]
  }
  # add %MT
  seur.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(seur.list[[i]], 
                                                          features = c("TrnP","TrnT","CYTB","TrnE","ND6","ND5","TrnL2",   
                                                                       "TrnS2","TrnH","ND4","ND4L","TrnR","ND3","TrnG",
                                                                       "COX3","ATP6","ATP8","TrnK","COX2","TrnD","TrnS1",
                                                                       "COX1","TrnY","TrnC","TrnN","TrnA","TrnW","ND2","TrnM",
                                                                       "TrnQ","TrnI","ND1","TrnL1","mt-Rnr2","TrnV","mt-Rnr1",
                                                                       "TrnF")) # Note: This is a list from ChrM from the REf Genome (mm39)
  
}


# See if any cells were removed:
cat((sum(unlist(lapply(mat.list, ncol)))-sum(unlist(lapply(seur.list, ncol)))),"cells (total) removed...\n")

## Preprocess seurat objects ##

## Function for following Seurat Workflow ##

seuPreProcess <- function(seu, assay='RNA', n.pcs=50, res=0.8){
  # NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')
  umap.name = paste0('umap_', assay)
  
  seu = NormalizeData(
    seu
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = F
  ) %>% ScaleData(
    assay = assay
  ) %>% RunPCA(
    assay = assay,
    reduction.name = pca.name,
    reduction.key = pca.key,
    verbose = F,
    npcs = n.pcs
  )
  
  # Find number of PCs #
  
  tmp.var <- (seu@reductions[[pca.name]]@stdev)^2
  var.cut <- 0.95*sum(tmp.var)
  j=0
  var.sum = 0
  while(var.sum < 0.95*var.cut){
    j = j + 1
    var.sum <- var.sum + tmp.var[j]
  }
  n.pcs.use = j    # Number of PCs to use for FindNeighbors 
  
  # FindNeighbors %>% RunUMAP, FindClusters
  
  seu <- FindNeighbors(
    seu,
    reduction = pca.name,
    dims = 1:n.pcs.use,
    force.recalc = TRUE,
    verbose = FALSE
  ) %>% RunUMAP(
    reduction = pca.name,
    dims = 1:n.pcs.use,
    reduction.name=umap.name
  )
  
  seu@reductions[[umap.name]]@misc$n.pcs.used <- n.pcs.use
  
  seu <- FindClusters(object = seu,resolution = res)
  seu[[paste0('RNA_res.',res)]] <- as.numeric(seu@active.ident)
  
  return(seu)
}  

## Preprocess Each Dataset Individually with Function ##

seur.list <- lapply(seur.list, seuPreProcess)

#### DoubletFinder on RNA for Eack Individual Dataset ####

bcmvn <- list()
pK <- list()
homotypic.prop <- list()
nExp_poi <- list()
nExp_poi.adj <- list()

## Prepare Doublet Estimation ##

for(i in 1:length(seur.list)){
  cat(' --------------------------------------------\n',
      '--- DoubletFinder for dataset number ', i, '---\n',
      '--------------------------------------------\n')
  
  ## pK Identification (no ground-truth)
  bcmvn[[i]]<- paramSweep_v3(             
    # pN-pK parameter sweeps on a 10,000-cell subset
    seu=seur.list[[i]], 
    PCs = 1:seur.list[[i]]@reductions$umap_RNA@misc$n.pcs.used, 
  ) %>% summarizeSweep(
    # computing the bimodality coefficient across pN and pK parameter space
    GT = FALSE
  ) %>% find.pK() 
  # Computes and visualizes the mean-variance normalized bimodality coefficient (BCmvn) score for each pK value tested during doubletFinder_ParamSweep
  
  ## Pull out max of bcmvn
  pK[[i]] <- as.numeric(as.character(bcmvn[[i]]$pK[bcmvn[[i]]$BCmetric==max(bcmvn[[i]]$BCmetric)])) # ugly, but functional...
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop[[i]] <- modelHomotypic(seur.list[[i]]$seurat_clusters) 
  
  nExp_poi[[i]] <- round(meta$`Multiplet Rate`[i]*length(colnames(seur.list[[i]])))  # Rate is 0.0242 for 5k "Target for Cell Recovery"
  nExp_poi.adj[[i]] <- round(nExp_poi[[i]]*(1-homotypic.prop[[i]]))
}

## Run DoubletFinder ##

for(i in 1:length(seur.list)){
  seur.list[[i]] <- 
    doubletFinder_v3( # just changed it so the output metadata column name is customizable
      seu=seur.list[[i]], 
      PCs = 1:seur.list[[i]]@reductions$umap_RNA@misc$n.pcs.used,
      pN = 0.25, #default value
      pK= pK[[i]], 
      nExp = nExp_poi.adj[[i]],  
      reuse.pANN = F
    )
}

#### Merging Datasets ####

# Add a prefix to the cell IDs so that cell IDs are NOT merged when the datasets are merged

Sample.ID <- c("mP2", "mP3", "mP4")

## Merge the datasets together ##

All.merged <- merge(seur.list[[1]], y = seur.list[2:length(Sample.ID)], add.cell.ids = Sample.ID[1:length(Sample.ID)])

## Save as RDS file ##

saveRDS(All.merged, file = "20220907_Mouse_Proximal.rds")

## Load RDS ##

All.merged <- readRDS(file = "20220907_Mouse_Proximal.rds", refhook = NULL)

#### Harmony Integration ####

# Need to normalize and run PCA before running Harmony with merged dataset

All.merged <- NormalizeData(object = All.merged, assay = 'RNA')
All.merged <- FindVariableFeatures(object = All.merged, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
All.merged <- ScaleData(object = All.merged, assay = 'RNA')
All.merged <- RunPCA(object = All.merged, assay = 'RNA', features = VariableFeatures(object = All.merged),
                     reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')

# Determine the 'dimensionality' of the dataset
ElbowPlot(All.merged,
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

## Define number of PCs ##

pcs <- npcs(seu = All.merged, var.toal = 0.95, reduction = 'pca_RNA')

## Run Harmony ##

gc() # clear some data

All.merged <- RunHarmony(All.merged, 
                         group.by.vars = c("Sample"), 
                         reduction = "pca_RNA", 
                         assay = "RNA", 
                         plot_convergence = TRUE)

## Normal Seurat Clustering ##

All.merged <- FindNeighbors(object = All.merged, 
                            reduction = "harmony",
                            k.param = 70, 
                            dims = 1:pcs, 
                            graph.name = 'harmony_snn')

All.merged <- FindClusters(object = All.merged, 
                           resolution = 0.9, 
                           graph.name = 'harmony_snn')


## Set Python UMAP via reticulate ##

umap.method = 'umap-learn'
metric = 'correlation'

All.merged <- RunUMAP(object = All.merged, reduction = "harmony", dims = 1:pcs)

## UMAP ##

DimPlot(object = All.merged, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

#### Visualizations for Proper Filtering ####

## Percent Mitochondria ##

VlnPlot(All.merged, features = c("percent.mt"), group.by = 'Sample', pt.size = 0,
        cols = natparks.pals(name="Glacier",n=3))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+     # Change X-Axis Text Size
  theme(axis.text.y = element_text(size = 16))+     # Change Y-Axis Text Size
  theme(axis.title.y = element_text(size = 18))+    # Change Y-Axis Title Text Size
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank())

## nFeature RNA ##

VlnPlot(All.merged, features = c("nFeature_RNA"), group.by = 'Sample', pt.size = 0,
        cols = natparks.pals(name="Glacier",n=3))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 18))+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank()) # Change object to visualize other samples

## nCount RNA ##

VlnPlot(All.merged, features = c("nCount_RNA"), group.by = 'Sample', pt.size = 0,
        cols = natparks.pals(name="Glacier",n=3))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 18))+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank()) # Change object to visualize other samples

## UMAP ##

DimPlot(object = All.merged, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = All.merged, 
        reduction = 'umap', 
        group.by = "Sample",
        cols = c('#81CAD6', '#EDCD44', '#DC3E26'),
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.3, 
        label.size = 5,
        shuffle = TRUE) + 
  labs(title = 'Colored by Sample')

## Feature Plot for Percent MT ##

FeaturePlot(All.merged, 
            reduction = "umap", 
            features = "percent.mt",
            pt.size = 0.4)+
  scale_color_viridis_c(option="F",begin=.2,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

## Feature Plot for nFeature RNA ##

FeaturePlot(All.merged, 
            reduction = "umap", 
            features = "nFeature_RNA",
            pt.size = 0.4)+
  scale_color_viridis_c(option="F",begin=.2,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

## Feature Plot for nCount RNA ##

FeaturePlot(All.merged, 
            reduction = "umap", 
            features = "nCount_RNA",
            pt.size = 0.4)+
  scale_color_viridis_c(option="F",begin=.2,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

## Feature Plot for specific Gene(s) ##

FeaturePlot(All.merged, 
            reduction = "umap", 
            features = "Epcam",
            pt.size = 0.2)+
  scale_color_viridis_c(option="F",begin=.2,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

#### Remove Low_Quality Cells ####

## Remove Cells w/ more than 30% MT genes / Less than 200 nFeature RNA / Less than 750 nCount RNA
# Note: These can come from literature and may vary based on tissue

Filtered_Cells <- subset(All.merged, subset = percent.mt < 30 &  
                           nFeature_RNA > 200 & nCount_RNA > 750 )

## Save as RDS file ##

saveRDS(Filtered_Cells, file = "20220907_Filtered_Cells.rds")

## Load RDS ##

Filtered_Cells <- readRDS(file = "20220907_Filtered_Cells.rds", refhook = NULL)

#### Re-Run Seurat Workflow + Harmony ####

## Check UMAP to See Where Cells Were Removed ##

DimPlot(object = Filtered_Cells, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

## Normalize and Run PCA before Running Harmony ##

Filtered_Cells <- NormalizeData(object = Filtered_Cells, assay = 'RNA')
Filtered_Cells <- FindVariableFeatures(object = Filtered_Cells, assay = 'RNA',
                                       selection.method = 'vst', nfeatures = 2000)
Filtered_Cells <- ScaleData(object = Filtered_Cells, assay = 'RNA')
Filtered_Cells <- RunPCA(object = Filtered_Cells, assay = 'RNA', 
                         features = VariableFeatures(object = Filtered_Cells),
                         reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')

## Determine the 'dimensionality' of the dataset ##

ElbowPlot(Filtered_Cells,
          reduction = 'pca_RNA',
          ndims = 50) 

##  Calculate the number of PCs that contain some proportion (95%) of the variance ##

pcs <- npcs(seu = Filtered_Cells, var.toal = 0.95, reduction = 'pca_RNA')

## Run Harmony ##

gc()

Filtered_Cells <- RunHarmony(Filtered_Cells,
                             group.by.vars = c("Sample"),
                             reduction = "pca_RNA", 
                             assay = "RNA", 
                             plot_convergence = TRUE)

saveRDS(Filtered_Cells, file = "20220907_Filtered_Cells_Harmony.rds")

## Normal Seurat Clustering ##

Filtered_Cells <- readRDS(file = "20220907_Filtered_Cells_Harmony.rds", refhook = NULL)


Filtered_Cells <- FindNeighbors(object = Filtered_Cells, 
                                reduction = "harmony", 
                                k.param = 40, 
                                dims = 1:pcs, 
                                graph.name = 'harmony_snn')

Filtered_Cells <- FindClusters(object = Filtered_Cells, 
                               resolution = 0.7, graph.name = 'harmony_snn')

## Set Python UMAP via reticulate ##

umap.method = 'umap-learn'
metric = 'correlation'


Filtered_Cells <- RunUMAP(object = Filtered_Cells, 
                          reduction = "harmony", dims = 1:pcs)

## UMAP ##

DimPlot(object = Filtered_Cells, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

FeaturePlot(Filtered_Cells, features = "Epcam")

#### Remove Doublets ####

DF.name = colnames(Filtered_Cells@meta.data)[grepl("DF.classification", colnames(Filtered_Cells@meta.data))]

## Find a Better Function to Combine All DF.classifications in the Future ##

mycol = coalesce(Filtered_Cells@meta.data$DF.classifications_0.25_0.28_26 , 
                 Filtered_Cells@meta.data$DF.classifications_0.25_0.3_6 ,
                 Filtered_Cells@meta.data$DF.classifications_0.25_0.3_25)

Filtered_Cells@meta.data$Doublet = mycol

# Doublet Dataframe #
pANN <- coalesce(Filtered_Cells@meta.data$pANN_0.25_0.28_26 , 
                 Filtered_Cells@meta.data$pANN_0.25_0.3_6 ,
                 Filtered_Cells@meta.data$pANN_0.25_0.3_25)

doublet <- Filtered_Cells@meta.data$Doublet

cluster <-  Filtered_Cells@meta.data$harmony_snn_res.0.7

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

write.csv(Avg_pANN , "Avg_pANN.csv")
write.csv(Pct_doublet , "Pct_doublet.csv")

## Visualize Doublets ##

DimPlot(object = Filtered_Cells, 
        reduction = 'umap', 
        group.by = "Doublet",
        cols = c( "#ffb6c1", "#380b11"),
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.4, 
        order = c("Doublet","Singlet"),
        label.size = 5) + 
  labs(title = 'Detecting Doublets')

## Compare to UMAP ##

DimPlot(object = Filtered_Cells, 
        reduction = 'umap', 
        group.by = "seurat_clusters",
        #cols = c( "#ffb6c1", "#380b11"),
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.5, 
        # = c("Doublet","Singlet"),
        label.size = 5) + 
  labs(title = 'Colored by Cluster')

## Stop here if not removing doublets ##

saveRDS(Filtered_Cells, file = "20220907_Filtered_Cells_Final.rds") 

#### Additional Doublet Removal if Needed ####
## Extra Plots ##

FeaturePlot(Filtered_Cells, 
            reduction = "umap", 
            features = "Krt5",
            pt.size = 0.4)+
  scale_color_viridis_c(option="F",begin=.2,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

DimPlot(object = Filtered_Cells, 
        reduction = 'umap', 
        group.by = "Location",
        cols = c('#ff99c8','#90e0ef','#9658B1'),
        repel = TRUE, 
        label = FALSE, 
        pt.size = 0.3, 
        label.size = 5,
        shuffle = TRUE) + 
  labs(title = 'Colored by Location')

## Remove Doublets ##

Idents(Filtered_Cells) <- 'sample' 

HQ_cells <- subset(
  Filtered_Cells,
  subset = Doublet == 'Singlet' 
)  


#### Re-Run Seurat Workflow + Harmony (w/out doublets) ####

## Need to Normalize and Run PCA before Running Harmony ##

HQ_cells <- NormalizeData(object = HQ_cells, assay = 'RNA')
HQ_cells <- FindVariableFeatures(object = HQ_cells, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
HQ_cells <- ScaleData(object = HQ_cells, assay = 'RNA')
HQ_cells <- RunPCA(object = HQ_cells, assay = 'RNA', features = VariableFeatures(object = HQ_cells),
                   reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')

## Determine the 'dimensionality' of the dataset ##

ElbowPlot(HQ_cells,
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

pcs <- npcs(seu = HQ_cells, var.toal = 0.95, reduction = 'pca_RNA')

## Run Harmony ##

gc()

HQ_cells <- RunHarmony(HQ_cells, 
                       group.by.vars = c("Patient","Process"), 
                       reduction = "pca_RNA", 
                       assay = "RNA", 
                       plot_convergence = TRUE)

## Normal Seurat clustering ##

HQ_cells <- FindNeighbors(object = HQ_cells, 
                          reduction = "harmony",
                          k.param = 70, 
                          dims = 1:pcs, 
                          graph.name = 'harmony_snn')

HQ_cells <- FindClusters(object = HQ_cells, 
                         resolution = 0.7, 
                         graph.name = 'harmony_snn')


## Set Python UMAP via reticulate ##

umap.method = 'umap-learn'
metric = 'correlation'

HQ_cells <- RunUMAP(object = HQ_cells, reduction = "harmony", dims = 1:pcs)

## UMAP ##

DimPlot(object = HQ_cells, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5) + 
  labs(title = 'Colored by Cluster')


## Save as RDS file ##

saveRDS(HQ_cells, file = "20220705_HQ_Cells.rds")


#### High Quality Cell Basic Analysis ####

## Load RDS ##

HQ_cells <- readRDS(file = "20220705_HQ_Cells.rds", refhook = NULL)

FeaturePlot(HQ_cells, features = c("KRT8"), pt.size = 0.75)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)
