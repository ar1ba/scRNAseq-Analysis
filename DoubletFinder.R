## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(filteredH486DFCm.seurat, PCs = 1:20, sct = FALSE)
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
annotations <- filteredH486DFCm.seurat$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(filteredH486DFCm.seurat@meta.data))  ## Assuming 7.6% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
filteredH486DFCm.seurat <- doubletFinder_v3(filteredH486DFCm.seurat, 
                                         PCs = 1:20, 
                                         pN = 0.25, #does not affect
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
#DimPlot(filteredH486DFCm.seurat, reduction = 'umap', group.by = "DF.classifications_0.25_0.28_1265")
#0.25 refers to pN, 0.28 is the pK value, and 1265 is nExp_poi.adj

# print counts of singlets and doublets
#table(filteredH486DFCm.seurat@meta.data$DF.classifications_0.25_0.28_1242)


#Filter out doublets
H486DFCm_clean <- subset(filteredH486DFCm.seurat, 
                           cells=rownames(filteredH486DFCm.seurat@meta.data)
                           [which(filteredH486DFCm.seurat@meta.data$DF.classification == "Singlet")])

#preprocess again

H486DFCm_clean <- NormalizeData(object = H486DFCm_clean)
H486DFCm_clean <- FindVariableFeatures(object = H486DFCm_clean)
H486DFCm_clean <- ScaleData(object = H486DFCm_clean)
H486DFCm_clean <- RunPCA(object = H486DFCm_clean) #linear dimensional data reduction

H486DFCm_clean <- FindNeighbors(H486DFCm_clean, dims = 1:20)

#choose cluster resolution
H486DFCm_clean <- FindClusters(H486DFCm_clean, resolution = c(0.04, 0.05, 0.1))
DimPlot(H486DFCm_clean, group.by = "RNA_snn_res.0.04", label = TRUE)
