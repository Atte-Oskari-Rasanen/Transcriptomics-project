ysis using Monocle
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

#run the script vai terminal by entering the script name (after you've saved it there)
# and add the rds object's name afterwards
#Pseudotime analysis using Monocle
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)
pdf(file="Monocl_pseudotime_subset_DI_FIN.pdf")
par(mfrow=(c(1,3)))
args = commandArgs(trailingOnly=TRUE)
rds_file <- args[1]

#source("order_cells.R")
singles <- readRDS(rds_file)

VlnPlot(singles, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 1, group.by = "assignment")

#subsetting applied. Different clusters taken in different subsets based on what batch correciton
#method was used. Initially the script was run without subsetting to see what the output looks lie
#after which the subsets were applied and the script was re-run

#For Harmony 
#Rigtclump <- singles[,  singles$seurat_clusters %in% c(17,13,12,4,5,1,3,2)]
#Leftclump <- singles[,  singles$seurat_clusters %in% c(0,18,23,21,10,7,19,17)] #maybe remove 25 and 14


#for DI
Rigtclump <- singles[,  singles$seurat_clusters %in% c(16,14,5,3,8,7)]
Leftclump <- singles[,  singles$seurat_clusters %in% c(20,15,18,17,1,21,9)]

#Assign each cluster to the most common cell type based on the made annotations
for(i in levels(singles)) {
  print(i)
  cells_to_reid <- WhichCells(singles, idents = i)
  newid <- names(sort(table(singles$seurat_clusters[cells_to_reid]),decreasing=TRUE))[1]
  Idents(singles, cells = cells_to_reid) <- newid
}




#define roots (PPC=pluripotent cell, since they should be treated as roots)
PPC <- singles[, singles$seurat_clusters == "14"]
head(PPC)
print("names")
PPC <- PPC$seurat_clusters
PPC<-names(PPC)


PPC2 <- singles[, singles$seurat_clusters == "15"]  # and 14 for DI
head(PPC2)
print("names")
PPC2 <- PPC2$seurat_clusters
PPC2<-names(PPC2)

#PPC
#assign clusters to cell types:


DefaultAssay(singles) <- "RNA"

##### plot without subsetting, convert seurat object into cell data set object required for monocle
singles.cds <- as.cell_data_set(singles)
singles.cds <- cluster_cells(cds = singles.cds, reduction_method = "UMAP")
singles.cds <- learn_graph(singles.cds, use_partition = TRUE)

singles1.cds <- order_cells(singles.cds, reduction_method = "UMAP", root_cells=PPC)
#singles2.cds <- order_cells(singles.cds, reduction_method = "UMAP", root_cells=PPC2)

DimPlot(singles, label = TRUE)




##### subsetting

Leftclump.cds <- as.cell_data_set(Leftclump)
Leftclump.cds <- cluster_cells(cds = Leftclump.cds, reduction_method = "UMAP")
Leftclump.cds <- learn_graph(Leftclump.cds, use_partition = TRUE)

Leftclump.cds <- order_cells(Leftclump.cds, reduction_method = "UMAP", root_pr_nodes=NULL,
                        root_cells=PPC2, verbose = FALSE)
#normally ppc2
print("left")
###plotting 
DimPlot(Leftclump, label = TRUE)

# plot trajectories colored by pseudotime
plot_cells(
  cds = Leftclump.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)



######
Rigtclump.cds <- as.cell_data_set(Rigtclump)
Rigtclump.cds <- cluster_cells(cds = Rigtclump.cds, reduction_method = "UMAP")
Rigtclump.cds <- learn_graph(Rigtclump.cds, use_partition = TRUE)

#print("right")
Rigtclump.cds <- order_cells(Rigtclump.cds, reduction_method = "UMAP", root_pr_nodes=NULL,
                           root_cells=PPC, verbose = FALSE)


####plotting
DimPlot(Rigtclump, label = TRUE)

# plot trajectories colored by pseudotime
plot_cells(
  cds = Rigtclump.cds,
  color_cells_by = "pseudotime",
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)




dev.off()


```
######################




#Pseudotime analysis using Slingshot. This was performed just for experimentation due to time 
#limit, did not produce meaningful results
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(SeuratWrappers)
args = commandArgs(trailingOnly=TRUE)
rds_file <- args[1]


singles <- readRDS("singles.combined.DI.rds")
pdf(file="Slingshot_DI_f2.pdf")
par(mfrow=(c(1,1)))


cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

VlnPlot(singles, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 1, group.by = "seurat_clusters")
DimPlot(singles, reduction = "pca",
        group.by = "seurat_clusters", pt.size = 2, label = TRUE, repel = TRUE) +
  scale_color_brewer(type = "qual", palette = "Set2")
DimPlot(singles, reduction = "umap",
        group.by = "seurat_clusters", pt.size = 2, label = TRUE, repel = TRUE) +
  scale_color_brewer(type = "qual", palette = "Set2")

sds <- slingshot(Embeddings(singles, "umap"), clusterLabels = singles$seurat_clusters,
                 start.clus = 15, stretch = 0)
print("cell pall")
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
print("cell col")
cell_colors <- cell_pal(singles$seurat_clusters, brewer_pal("qual", "Set2"))

cell_colors_clust <- cell_pal(singles$seurat_clusters, hue_pal())

#plot(reducedDim(singles), col = cell_colors, pch = 16, cex = 0.5)
#lines(singles, lwd = 2, type = 'lineages', col = 'black')

#plot(reducedDim(singles), col = cell_colors_clust, pch = 16, cex = 0.5)
#lines(singles, lwd = 2, type = 'lineages', col = 'black')


#add colours to the clusters
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.7)
lines(sds, lwd = 2, type = 'lineages', col = 'black')

#principal curves version
plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.7)
lines(sds, lwd = 2, col = 'black')


####which cell is in which lineage?
nc <- 3
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(50, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.7, main = i)
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}
dev.off()


```
######################
