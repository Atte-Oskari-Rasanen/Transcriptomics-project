#3a. Batch effect correction using Harmony
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

###Harmony step (comes after step 2)
pdf(file='harmony_2015.pdf')
par(mfrow=(c(1,3)))

library(cowplot)
library(harmony)
library(Rcpp)
#the scripts below had to be sourced to get them working
source('do_scatter.R')
source('RunHarmony.R')
#modify clustering dims and res as you see fit
singles <- RunHarmony(singles, 'assignment')
singles <- FindNeighbors(singles, dims = 1:25)#builds the shared nearest neighbourhood (SNN) graph for the dataset
singles <- FindClusters(singles, resolution = 1.5) #cluster identification based on the SNN
singles <- RunUMAP(singles, reduction = "harmony", dims = 1:25) 

#plotting
DimPlot(singles, reduction = "harmony", label = T)
FeaturePlot(singles, features = 'GFP') -> plot1
FeaturePlot(singles, features = 'EBFP') -> plot2
FeaturePlot(singles, features = 'CHERRY') -> plot3
DimPlot(singles, group.by = 'assignment') -> plot4 

plot_grid(plot1, plot2, plot3, plot4, ncol = 2)

dev.off()
### Finding differentially expressed features
#find markers, apply thresholds to filter out poor genes
singles.unstranded.markers <- FindAllMarkers(singles, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = 0.4)
t <- table(Idents(singles), singles$assignment)

#write into an output file the cluster x assignment table for later
write.csv(t, '/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratOutputs/20_15/cluster_x_assignment_20_15.csv')


#write into an output file the marker genes
#into csv
write.csv(singles.unstranded.markers, '/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratOutputs/20_15/markers.harmony_1515.csv')




```
######################




#3a.2 If you want to apply harmony function rather in a loop:
######################
```{r setup, include=FALSE}
args = commandArgs(trailingOnly=TRUE)
#here looping over harmony and trying out different options of clustering parameters

disMax <- 40 
resMax <- 3
dims <- dimsMax[seq(10, 40, 5)] #so going over every 5th element starting from 10 and ending in 40
res <- resMax[seq(0.5,3,0.5)] #same logic 

source('do_scatter.R')
source('RunHarmony.R')

#Function for batch effect correction and generation of umap plots:

harmony_and_umap <- function(dim, res){
###Harmony step
pdf(file='harmony_', d, '_', r, '.pdf')
par(mfrow=(c(1,3)))

library(cowplot)
library(harmony)
library(Rcpp)
source('do_scatter.R')
source('RunHarmony.R')

singles <- RunHarmony(singles, 'assignment')
singles <- FindNeighbors(singles, dims = 1:dim)
singles <- FindClusters(singles, resolution = res)
singles <- RunUMAP(singles, reduction = "harmony", dims = 1:dim) 

DimPlot(singles, reduction = "harmony", label = T)
FeaturePlot(singles, features = 'GFP') -> plot1
FeaturePlot(singles, features = 'EBFP') -> plot2
FeaturePlot(singles, features = 'CHERRY') -> plot3
DimPlot(singles, group.by = 'assignment') -> plot4 

plot_grid(plot1, plot2, plot3, plot4, ncol = 2)

dev.off()

### Finding differentially expressed features
#find markers
singles.unstranded.markers <- FindAllMarkers(singles, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = 0.4)
t <- table(Idents(singles), singles$assignment)

#t <- table(Idents(singles), singles$assignment)
write.csv(t, '/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratOutputs/20_15/cluster_x_assignment_', dim, '_', res, '.csv')


#write into an output file
#into csv
write.csv(singles.unstranded.markers, '/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratOutputs/20_15/markers.harmony_', dim, '_', res, '.csv')

######

singles <- FindNeighbors(singles, dims = 1:dim)
singles <- FindClusters(singles, resolution = res)

head(Idents(singles), 5)

singles <- RunUMAP(singles, dims = 1:dim)

pdf(file="plots3_allclusters_', dim, '_', res, '.pdf")
par(mfrow=(c(1,3)))


DimPlot(singles, reduction = "umap", label = T)
FeaturePlot(singles, features = 'GFP') -> plot1
FeaturePlot(singles, features = 'EBFP') -> plot2
FeaturePlot(singles, features = 'CHERRY') -> plot3
DimPlot(singles, group.by = 'assignment') -> plot4 

library(cowplot)

#plot all 4 side by side for viewing
plot_grid(plot1, plot2, plot3, plot4, ncol = 2)
dev.off()
}


for (d in dims){
  for(r in res){
    harmonyfunction(d,r)
  }
} 
```
######################




#3b. Batch effect correction using Seurat (comes after step 2)
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

##### Data integration
library(ggplot2)
#DimPlot(singles, label = T)
#now try integration
# split the dataset into a list of three seurat objects 
singles.list <- SplitObject(singles, split.by = "assignment")

# normalize and identify variable features for each dataset independently
#changed to 8000 variable features
singles.list <- lapply(X = singles.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 8000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = singles.list)

# make the anchors based on the features
singles.anchors <- FindIntegrationAnchors(object.list = singles.list, anchor.features = features)


###Creating an integrated data assay
# this command creates an 'integrated' data assay
singles.combined <- IntegrateData(anchorset = singles.anchors)

print(singles.combined)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(singles.combined) <- "integrated"

# Run the standard workflow for visualization and clustering

pdf('DI_plots_25_15_umap.pdf')

par(mfrow=(c(1,3)))

# Run the standard workflow for visualization and clustering
singles.combined <- ScaleData(singles.combined, verbose = FALSE)
singles.combined <- RunPCA(singles.combined, npcs = 25, verbose = FALSE)
singles.combined <- RunUMAP(singles.combined, reduction = "umap", dims = 1:25)
singles.combined <- FindNeighbors(singles.combined, reduction = "umap", dims = 1:25)
singles.combined <- FindClusters(singles.combined, resolution = 1.5)


DimPlot(singles.combined, label = T)
DimPlot(singles.combined, group.by = 'assignment')
DimPlot(singles.combined, group.by = 'assignment') + facet_wrap(~singles.combined$assignment, ncol=3)

#do a little QC
VlnPlot(singles.combined, features = 'percent.mito') 
VlnPlot(singles.combined, features = 'nFeature_RNA') 
VlnPlot(singles.combined, features = 'nCount_RNA')
dev.off()

#write outputs same as with harmony
singles.combined.markers <- FindAllMarkers(singles.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = 0.4)
#write.csv(singles.combined.markers, paste0('DI_singles_integr_markers',d,'_',r,'.csv'))
t <- table(Idents(singles.combined), singles.combined$assignment)
#write.csv(t, paste0('DI_cluster_x_assignment_',d,'_',r,'.csv'))


```
######################




#Earlier, findmarkers function was used but decided to run findconservedmarkers separately as well
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

library(Seurat)
library(tximport)
library(dplyr)
library(devtools)
#script for applying find conserved markers for each output file
#as an input r object need to be either the one where data integration (seurat batch correction)
#via seurat has been performed or the one with harmony.
#the settings change based on the number of dimensions used

#upload the object where the data was batch corrected (in this example via harmony)
singles.combined <- readRDS("singles.harmony.rds")

DefaultAssay(singles.combined) <- "RNA"

#ident.1 stands for the cluster for which the markers should be defined for 
for (i in 1:23){
#Idents(x = singles.combined)
skip_to_next <- FALSE

#apply findocnservedmarkers function, put it inside trycatch to skip error messages etc.
# the error messages arise when there is a file where no conserved markers could be found
tryCatch(singles.consv <- FindConservedMarkers(
  singles.combined,
  ident.1 = i,
  ident.2 = NULL,
  grouping.var = "assignment",
  assay = "RNA",
  slot = "data",
  meta.method = metap::minimump,
  verbose = TRUE,
)
, error = function(singles.consv) { skip_to_next <<- TRUE}) #skip an error without crashing the script

if(skip_to_next) { next }
saveRDS(singles.combined, file='noBESingles.with.clustr2515.rds')
write.csv(singles.consv, paste0('singles.conserv.h',i,'.csv')) #output the cluster data
}


```
######################


