
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

singles.anchors <- FindIntegrationAnchors(object.list = singles.list, anchor.features = features)


###Creating an integrated data assay
# this command creates an 'integrated' data assay
singles.combined <- IntegrateData(anchorset = singles.anchors)

print(singles.combined)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(singles.combined) <- "integrated"

# Run the standard workflow for visualization and clustering

# Run the standard workflow for visualization and clustering
singles.combined <- ScaleData(singles.combined, verbose = FALSE)
singles.combined <- RunPCA(singles.combined, npcs = 25, verbose = FALSE)
singles.combined <- RunUMAP(singles.combined, reduction = "pca", dims = 1:25)
singles.combined <- FindNeighbors(singles.combined, reduction = "pca", dims = 1:25)
singles.combined <- FindClusters(singles.combined, resolution = 1.5)

#dev.off()
# pdf(paste0(file='DI_plots_', d, '_', r, '.pdf'))
# par(mfrow=(c(2,3)))
# 
 DimPlot(singles.combined, label = T)
 DimPlot(singles.combined, group.by = 'assignment')
 DimPlot(singles.combined, group.by = 'assignment') + facet_wrap(~singles.combined$assignment, ncol=3)
# 
singles.combined <- RenameIdents(object = singles.combined, `0` = "0", `1` = "Developing T cell", `2` = "2", `3`= "immature/prolif B cell",
                        `4`="4", `5`="B cell", `6`="6", `7`="Immature macrophage", `8`= "immture Macrophage",
                        `9`="Prolif macrophage", `10`="10", `11`= "11", `12`="B cell",`13`="13", `14`="Proliferating HPC",
                        `15`="Immature mast cell/differentiating HPC", `16`="Erythrocyte", `17`="Erythrocyte", `18`="Mast cell", 
                        `20`="prolif. B cell")


 
 saveRDS(singles.combined, file = 'singles.combined.DI.rds')

 dev.off()
```



