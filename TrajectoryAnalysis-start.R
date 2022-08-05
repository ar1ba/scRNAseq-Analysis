##INSTALL PACKAGES
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)

# ...1 Convert to cell_data_set object ------------------------

library(Seurat)
library(SeuratDisk)
library(SeuratData)
allsamples <- LoadH5Seurat("~/Downloads/allsamples.h5seurat")

cds <- as.cell_data_set(allsamples)
cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = "louvain")
colData(cds)
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
counts(cds)


cds <- preprocess_cds(cds, method = "PCA")
cds = align_cds(cds, num_dim = 50, alignment_group = "plate")

cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "ident")

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# ...5. Finding genes that change as a function of pseudotime --------------------
deg_allcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_allcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

FeaturePlot(cds, features = c('HES4', 'ISG15', 'GNB1', "GABRD", "SKI", "HES5"))


# visualizing pseudotime in seurat

allsamples$pseudotime <- pseudotime(cds)
Idents(allsamples) <- allsamples$seurat_clusters
FeaturePlot(allsamples, features = "pseudotime", label = T)
