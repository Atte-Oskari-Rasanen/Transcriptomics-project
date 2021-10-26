```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

#the scripts were kept separate while working which is why sourcing the earlier components here,
#can also upload the object if the computer is slow:
#singles <- readRDS(singles.rds)
source("1_SeuratObjectCreation.R")

#remove outliers based on the earlier created scatter plots
singles <- subset(singles, subset = nFeature_RNA > 500 & nFeature_RNA < 12000 & percent.mito < 0.06)
#normalising the data via log transofrmation 
singles <- NormalizeData(singles, normalization.method = "LogNormalize", scale.factor = 10000)


##############
singles <- FindVariableFeatures(singles, selection.method = "vst", nfeatures = 8000)
#singles_mvp <- FindVariableFeatures(singles, selection.method = "mvp", nfeatures = 8000)

#singles_disp <- FindVariableFeatures(singles, selection.method = "disp", nfeatures = 8000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(singles), 10)

write.csv(top10, '/home/inf-54-2020/trans_aligned_alevin_15000_output_transindex/SeuratOutputs/20_15/top10variablegenes_2015.csv')

# plot variable features with and without labels
pdf(file="plots2.pdf")
par(mfrow=(c(1,3)))


plot1 <- VariableFeaturePlot(singles)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #label the top10 genes
plot1 + plot2

#scale the data of singles
singles <- ScaleData(singles)


#perform pca analysis based on variablefeatures 
singles <- RunPCA(singles, features = VariableFeatures(object = singles))

print(singles[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(singles, dims = 1:2, reduction = "pca")

DimPlot(singles, reduction = "pca")

DimHeatmap(singles, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(singles, dims = 1:25, cells = 500, balanced = TRUE)

###deciding on the numbner of dimensions for clustering
#Elbowplot

ElbowPlot(singles, ndims = 50)

#Trying out jackstraw as well for determining the appropriate number of PCs for clustering
pbmc <- JackStraw(singles, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pbmc
JackStrawPlot(pbmc, dims = 1:15)

dev.off()

```
######################



#Scripts that needed to be sourced to run Harmony. Normally would be part of Harmony package but had to
#source them separately due to unknown error messages
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

#do_scatter found from https://github.com/immunogenomics/harmony/blob/master/vignettes/quickstart.Rmd

do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE,
                       do_labels = TRUE, nice_names, 
                       palette_use = colors_use,
                       pt_size = 4, point_size = .5, base_size = 12, 
                       do_points = TRUE, do_density = FALSE, h = 6, w = 8) {
    umap_use <- umap_use[, 1:2]
    colnames(umap_use) <- c('X1', 'X2')
    plt_df <- umap_use %>% data.frame() %>% 
        merge(meta_data) %>% 
        dplyr::sample_frac(1L) 
    plt_df$given_name <- plt_df[[label_name]]
    
    if (!missing(nice_names)) {
        plt_df %<>%
            dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))
        
        plt_df[[label_name]] <- plt_df$nice_name        
    }
    
    plt <- plt_df %>% 
        ggplot(aes_string("X1", "X2", col = label_name, fill = label_name)) + 
        theme_test(base_size = base_size) + 
        theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1,
                                                        shape = 16, size = 4)), 
               alpha = FALSE) +
        scale_color_manual(values = palette_use) + 
        scale_fill_manual(values = palette_use) +    
        theme(plot.title = element_text(hjust = .5)) + 
        labs(x = "PC 1", y = "PC 2") 
    
    if (do_points) 
        plt <- plt + geom_point(shape = '.')
    if (do_density) 
        plt <- plt + geom_density_2d()    
    
    
    if (no_guides)
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
    
    if (do_labels) {
        data_labels <- plt_df %>% 
            dplyr::group_by_(label_name) %>% 
            dplyr::summarise(X1 = mean(X1), X2 = mean(X2)) %>% 
            dplyr::ungroup()
        
        plt <- plt + geom_label(data = data_labels, label.size = NA,
                        aes_string(label = label_name), 
                        color = "white", size = pt_size, alpha = 1,
                        segment.size = 0) +
                guides(col = FALSE, fill = FALSE)
    }
    
    return(plt)
}


