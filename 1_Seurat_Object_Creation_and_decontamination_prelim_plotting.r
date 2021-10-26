#Data processing in Seurat ( on R ver 4.0.5). The paramateres used are example ones and should be changed according to own analysis 
######################
#1. Seurat Object creation and decontamination
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

library(Seurat)
library(tximport)
library(dplyr)
library(devtools)

#set working dir. according to the one where the folder and file structure is (as described above)
setwd('~/trans_aligned_alevin_15000_output_transindex')
#load in data
files <- file.path('~/trans_aligned_alevin_15000_output_transindex/alevin/quants_mat.gz')
file.exists(files)

#importing all the files (json)
txi <- tximport(files, type="alevin")

#Focus on barcode, status, assignment
#load metadata from souporcell
meta <- read.table('~/trans_aligned_alevin_15000_output_transindex/alevin/modified.input.clusters.tsv', sep = '\t', header = T,row.names = 1)

#creating a seurat object
spleen.soup <- CreateSeuratObject(counts = txi$counts, meta.data = meta, min.cells = 3, min.features = 200, project = "10X_splenocytes_Alevin_IU")
spleen.soup



pdf(file="/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratPlots/plots1.pdf")
par(mfrow=(c(1,3)))

####Cleaning the dataset - convert the seurat object first to a singlecellexperiemnt one, 
#otherwise decontX will not work!
library(celda)
sce <- as.SingleCellExperiment(spleen.soup)
sce
#run decontx, let it cluster
sce <- decontX(x = sce)
sce$decontX_contamination
plotDecontXContamination(sce)
decontXcounts(sce) -> counts
#create a cleaned up seurat object
spleen.decon <- CreateSeuratObject(counts = counts, meta.data = meta, min.cells = 3, min.features = 200, project = "10X_splenocytes_Alevin_Ievin_IU")

saveRDS(spleen.decon, file = 'spleen.soup.decon.rds')

```
######################



#1.2 Scanning for mitochondrial genes and preliminary plotting of the data
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

spleen.decon <- readRDS('spleen.soup.decon.rds')
#continuing with subset singlets
subset(spleen.decon, subset = status == 'singlet') -> singles
#match_df(df1, 0, on = NULL)


#find mito genes
grep(rownames(singles), pattern = '-COX1-', value = T) -> Cox1
grep(rownames(singles), pattern = '-COX2-', value = T) -> Cox2
grep(rownames(singles), pattern = '-ATP8-', value = T) -> Atp8
grep(rownames(singles), pattern = '-ATP6-', value = T) -> Atp6
grep(rownames(singles), pattern = '-COX3-', value = T) -> Cox3
grep(rownames(singles), pattern = '-NU1M-', value = T) -> nu1m
grep(rownames(singles), pattern = '-NU2M-', value = T) -> nu2m
grep(rownames(singles), pattern = '-NU3M-', value = T) -> nu3m
grep(rownames(singles), pattern = '-NU4M-', value = T) -> nu4m
grep(rownames(singles), pattern = '-NU4LM-', value = T) -> nu4lm
grep(rownames(singles), pattern = '-NU5M-', value = T) -> nu5m
grep(rownames(singles), pattern = '-NU6M-', value = T) -> nu6m
grep(rownames(singles), pattern = '-CYB-', value = T) -> cyb

#add these all together - put into a list
c(Cox1,Cox2,Atp8,Atp6,Cox3,nu1m,nu2m,nu3m,nu4m,nu5m,nu6m,cyb) -> mito.genes

#make percent.mito
percent.mito <- Matrix::colSums(GetAssayData(singles, slot = "counts")[mito.genes,])/Matrix::colSums(GetAssayData(singles, slot = "counts"))

#colsums computes the sums of the matrix, takes singles
#getAssaData pulls info for specified stored dimensional reduction analysis, slot specifies
#the kind of data to retrieve, singles is the(seurat) object

#add to seurat object
singles$percent.mito <- percent.mito

#get a violin plot to check the data
VlnPlot(singles, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

#check the spread of the data
plot1 <- FeatureScatter(singles, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(singles, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2

saveRDS(singles, file = 'singles.rds') #save the object for later

dev.off()

```
######################


