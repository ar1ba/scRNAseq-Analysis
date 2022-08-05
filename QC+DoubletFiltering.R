
#load libraries ---------------------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)

# load matrix from shared folder ---------------------------------------------------------------------------------------

H488DFCm <- Read10X(data.dir = "/gpfs/ycga/project/sestan/sestanlabShare_from_Dan/CellRanger/H488DFCm/RNAalone/filtered_feature_bc_matrix")
H488DFCm.seurat.obj <- CreateSeuratObject(counts = H488DFCm, project = 'matrix', min.cells = 
                                          3, min.features = 200)

H488DFCm.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(H488DFCm.seurat.obj,pattern = "^MT-")

# create a filtered seurat object, separate from raw data; subset by mitochondrial percentage, nCount_RNA, and nFeatureRNA
# visualize with violin plot
VlnPlot(filteredH488DFCm.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

filteredH488DFCm.seurat <- subset(H488DFCm.seurat.obj, subset = percent.mt < 5 & nCount_RNA > 200 & nFeature_RNA > 200)

# standard preprocessing workflow ---------------------------------------------------------------------------------------
filteredH488DFCm.seurat <- NormalizeData(object = H488DFCm.seurat.obj)
filteredH488DFCm.seurat <- FindVariableFeatures(object = filteredH488DFCm.seurat)
filteredH488DFCm.seurat <- ScaleData(object = filteredH488DFCm.seurat)
filteredH488DFCm.seurat <- RunPCA(object = filteredH488DFCm.seurat) #linear dimensional data reduction

# visualize PCA's in elbow plot and determine dimensions ---------------------------------------------------------------------------------------
ElbowPlot(filteredH488DFCm.seurat)

filteredH488DFCm.seurat <- FindNeighbors(object = filteredH488DFCm.seurat,
                                       dims = 1:20 )
filteredH488DFCm.seurat <- FindClusters(object = filteredH488DFCm.seurat)

# plot clusters on a UMAP ---------------------------------------------------------------------------------------
filteredH488DFCm.seurat <- RunUMAP(object = filteredH488DFCm.seurat, dims = 1:20)


## Begin identifying doublets using DoubletFinder (R) ---------------------------------------------------------------------------------------
library(DoubletFinder)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(filteredH488DFCm.seurat, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn_matrix <- find.pK(sweep.stats)

ggplot(bcmvn_matrix, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

#save optimal pK value
pK <- bcmvn_matrix %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- filteredH488DFCm.seurat$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)          
nExp_poi <- round(0.076*nrow(filteredH488DFCm.seurat@meta.data))  ## Assuming 7.6% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder ---------------------------------------------------------------------------------------
filteredH488DFCm.seurat <- doubletFinder_v3(filteredH488DFCm.seurat, 
                                          PCs = 1:20, 
                                          pN = 0.25, #does not affect
                                          pK = pK, 
                                          nExp = nExp_poi.adj,
                                          reuse.pANN = FALSE, sct = FALSE)


# visualize doublets (adjust classificiations based on pK and nExp_poi.adj -------------------------------------------------------------
DimPlot(filteredH488DFCm.seurat, reduction = 'umap', group.by = "DF.classifications_0.25_0.24_710")

# list number of singlets and doublets
table(filteredH488DFCm.seurat@meta.data$DF.classifications_0.25_0.24_710)


#Filter out doublets and create clean Seurat object ---------------------------------------------------------------------------------------
H488DFCm_clean <- subset(filteredH488DFCm.seurat, 
                       cells=rownames(filteredH488DFCm.seurat@meta.data)
                       [which(filteredH488DFCm.seurat@meta.data$DF.classification == "Singlet")])
