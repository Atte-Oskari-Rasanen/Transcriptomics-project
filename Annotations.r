#Creating the annotated seurat objects. Annotation of the clusters was done manually based on
#available literature


#adding the annotations
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

library(Seurat)
library(tximport)
library(dplyr)
library(devtools)

singles <- readRDS('singles.rds')

library(cowplot)
library(harmony)
library(Rcpp)
library(ggplot2)

source('do_scatter.R')
source('RunHarmony.R')
#harmony_and_umap <- function(d, r){
###Harmony step
#pdf(paste0("harmony.pdf"))
#jpeg(paste0("harmony",d, "_", r ,".jpeg", sep="")
#par(mfrow=(c(1,3)))

singles <- RunHarmony(singles, 'assignment')
singles <- FindNeighbors(singles, dims = 1:25)
singles <- FindClusters(singles, resolution = 1.5)
singles <- RunUMAP(singles, reduction = "harmony", dims = 1:25) 


singles <- RenameIdents(object = singles, `0` = "0", `1` = "1", `2` = "2", `3`= "immature/prolif B cell", `4`="4",
                        `5`="T cell", `6`="Immature macrophage", `7`="7", `8`= "Mature Macrophage", `9`="9", 
                        `10`="SMC", `11`= "11", `12`="Immature/prolif T cell/NK",`13`="Differentiating cell", `14`="14",
                        `15`="Neutrophil", `16`="16", `17`="Erythrocyte", `18`="CD8 T cell", `19`="Immature Mast cell", 
                        `20`="20", `21`="Lymphocyte", `22`="Helper T cell", `23`="Immature Macrophage", `24`="24", `25`="25")

singles$seurat_clusters
saveRDS(singles, file = 'singles.harmony.rds')

```


#singles with seurat batch correction
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

library(Seurat)
library(SeuratData)
library(patchwork)

library(tximport)
library(dplyr)
library(devtools)
library(celda)
library(SeuratWrappers)

singles <- readRDS('singles.rds')


#colsums computes the sums of the matrix, takes singles
#getassadata pulls info for specified stored dimensional reduction analysis, slot specifies
#the kind of data to retrieve, singles is the(seurat) object

#add to seurat object
singles$percent.mito <- percent.mito
singles <- subset(singles, subset = nFeature_RNA > 500 & nFeature_RNA < 12000 & percent.mito < 0.06)


pdf(file="/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratPlots/DI_plots1_2515.pdf")
par(mfrow=(c(1,3)))
