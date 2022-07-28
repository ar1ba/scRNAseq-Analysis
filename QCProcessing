#Human Sample 486, Dorsolateral Frontal Cortex
#Sex: Male	PMI Hours: 22	pH: 6.21	RIN: 6	Age: 0.3	Years	
# AgeDays: 375.5	Period: 8	Epoch: 2

#load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

# Load in Matrix from 10x Cell Ranger Results
H486DFCm <- Read10X(data.dir = "/gpfs/ycga/project/sestan/sestanlabShare_from_Dan/CellRanger/H486DFCm/RNAalone/filtered_feature_bc_matrix")

#Create Seurat Object
H486DFCm.seurat.obj <- CreateSeuratObject(counts = H486DFCm, project = 'matrix', min.cells = 
                                          3, min.features = 200)

#QC
#Mitochondrial Percent filtering

H486DFCm.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(H486DFCm.seurat.obj,pattern = "^MT-")
filteredH486DFCm.seurat <- subset(H486DFCm.seurat.obj, subset = percent.mt < 5)

# standard preprocessing workflow
filteredH486DFCm.seurat <- NormalizeData(object = H486DFCm.seurat.obj)
filteredH486DFCm.seurat <- FindVariableFeatures(object = filteredH486DFCm.seurat)
filteredH486DFCm.seurat <- ScaleData(object = filteredH486DFCm.seurat)
filteredH486DFCm.seurat <- RunPCA(object = filteredH486DFCm.seurat) #linear dimensional data reduction
ElbowPlot(filteredH486DFCm.seurat)
filteredH486DFCm.seurat <- FindNeighbors(object = filteredH486DFCm.seurat,
                                       dims = 1:20 )
filteredH486DFCm.seurat <- FindClusters(object = filteredH486DFCm.seurat)
filteredH486DFCm.seurat <- RunUMAP(object = filteredH486DFCm.seurat, dims = 1:20)
