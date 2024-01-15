 ##################################################################
#  (1a) No regression:  #
##################################################################

rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(future)
plan("multiprocess", workers = 256)
options(future.globals.maxSize = 256000 * 1024^2)

#PDX.combined <- readRDS("/mnt/Sara/write/Danielli_Patel_Langenau_20221220.rds")
PDX.combined <- readRDS("/Volumes/Sara_PhD/scRNAseq_data/write/Danielli_Patel_Langenau_20221220.rds")
head(PDX.combined@meta.data)

## Data normalization
PDX.combined <- NormalizeData(PDX.combined, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable features
PDX.combined <- FindVariableFeatures(PDX.combined, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(PDX.combined), 10)
# Scaling data
PDX.combined <- ScaleData(PDX.combined, features = VariableFeatures(PDX.combined))

## Dimension reduction
PDX.combined
PDX.combined <- RunPCA(PDX.combined, reduction.name = "pca_no_regression", reduction.key = "pca_no_regression")

PDX.combined$origin <- factor(x = PDX.combined$origin, levels = c("Danielli et al.", "Patel et al.", "Wei et al.", "Cheng et al."))

PDX.combined$name <- factor(x =PDX.combined$name, levels = c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3', 
                                                             'aRMS-4',  'aRMS-5', 'KFR', 'Mast118', 
                                                             'Mast95', 'MSK72117', 'MSK72117_SC', 
                                                             'MSK82489', 'Rh4',  'Rh41',   'RMS', 
                                                             'SJRHB010468_D1',  'SJRHB010468_X1', 'SJRHB013757_D2', 
                                                             'SJRHB013757_X1', 'SJRHB013759_A1',  'SJRHB013759_A2', 
                                                             'SJRHB013759_X14','SJRHB013759_X15', 'SJRHB031320_D1', 
                                                             'SJRHB031320_X1', 'SJRHB046156_A1', 'SJRHB046156_X1', 
                                                             '20696','21202', '29806', 
                                                             'eRMS-1.1','eRMS-1.2',  'eRMS-2.1', 'eRMS-2.2', 
                                                             'eRMS-3.2','eRMS-4',  'eRMS-8.1', 'eRMS-8.2', 
                                                             'eRMS-8.3', 'Mast111','Mast139',  'Mast139_SC', 
                                                             'Mast39',   'Mast85_r1','Mast85_r2', 
                                                             'Mast85_r2_SC', 'MSK74711', 'RD', 'SJRHB000026_R2',  'SJRHB000026_R3',  'SJRHB000026_X1', 
                                                             'SJRHB000026_X2', 'SJRHB010927_D1', 'SJRHB010927_X1', 
                                                             'SJRHB010928_R1', 'SJRHB010928_X1',  'SJRHB011_D', 
                                                             'SJRHB011_X', 'SJRHB012_R', 'SJRHB012_S', 'SJRHB012_Y', 
                                                             'SJRHB012_Z', 'SJRHB012405_D1', 'SJRHB012405_X1',  'SJRHB013758_D1', 
                                                             'SJRHB013758_D2', 'SJRHB013758_X1', 'SJRHB013758_X2',  'SJRHB030680_R1', 
                                                             'SJRHB030680_X1', 'SJRHB049189_D1', 'SJRHB049189_X1'))


## Color scales
# Color model
paletteer::paletteer_d("rcartocolor::Tropic", n = 4)
col1 <- c("#009B9EFF","#A7D3D4FF",  "#E4C1D9FF","#C75DABFF")

# Color origin
col2 <- paletteer::paletteer_d("ggthemes::excel_Slice", n=4)

# Color name
col3 <- paletteer::paletteer_c("ggthemes::Blue-Teal", n = 27)
col4 <- paletteer::paletteer_c("ggthemes::Purple", n = 45)




DimPlot(PDX.combined, reduction = "pca_no_regression", group.by = 'name', cols = c(col3, col4))
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/1_PCA_name.pdf", 
       width=13, height=5, dpi=300)

DimPlot(PDX.combined, reduction = "pca_no_regression", group.by = 'model', cols = col1, order = c('Cell line', 'Primary culture', 'O-PDX', 'Patient'),
        pt.size = 0.1) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/1_PCA_model.pdf", 
       width=6, height=5, dpi=300)


ElbowPlot(PDX.combined, reduction='pca_no_regression')
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/2_Elbowplot.pdf", 
       width=6, height=4, dpi=300)


## Cluster cells
PDX.combined <- FindNeighbors(PDX.combined, reduction = 'pca_no_regression', dims = 1:15)

# UMAP
PDX.combined <- RunUMAP(PDX.combined, dims = 1:15, reduction = 'pca_no_regression', reduction.name = 'UMAP_no_regression')

## Plots
# Umap plot based on sample
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'name', pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE,
        cols = c(col3, col4))+ NoLegend() + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_col_noleg.pdf", 
       width=5, height=5, dpi=300)

DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'name', pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE,
        cols = c(col3, col4)) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_col.pdf", 
       width=15, height=5, dpi=300)

## Umap plot based on RMS subtype
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'subtype', pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE,
        order = 'aRMS', cols = c('#cc96b7ff', '#4aa2bcff')) + NoLegend() + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_subtype_noleg.pdf", 
       width=5, height=5, dpi=300)


## Umap plot based on RMS model
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'model',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE, cols = col1, order = c('Cell line', 'Primary culture', 'O-PDX', 'Patient')) + NoLegend() + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_model_noleg.pdf", 
       width=5, height=5, dpi=300)


## UMAP plot by origin
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'origin', order = 'Cheng et al',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE, cols = col2)  + NoLegend() + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_origin_noleg.pdf", 
       width=5, height=5, dpi=300)


## UMAP plot by site
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'site',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()+ NoLegend()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_site_noleg.pdf", 
       width=5, height=5, dpi=300)

## UMAP plot by status
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'status',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()+ NoLegend()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_status_noleg.pdf", 
       width=5, height=5, dpi=300)





### Plots with legend

# Umap plot based on sample
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'name', pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE,
        cols = c(col3, col4))
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_col.pdf", 
       width=13, height=5, dpi=300)


## Umap plot based on RMS subtype
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'subtype', pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE,
        order = 'aRMS', cols = c('#cc96b7ff', '#4aa2bcff')) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_subtype.pdf", 
       width=5, height=4, dpi=300)

## Umap plot based on RMS model
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'model',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE, cols = col1, order = c('Cell line', 'Primary culture', 'O-PDX', 'Patient'))  + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_model.pdf", 
       width=5, height=4, dpi=300)

## UMAP plot by origin
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'origin', order = 'Cheng et al',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE, cols = col2) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_origin.pdf", 
       width=5, height=4, dpi=300)

## UMAP plot by site
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'site',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_site.pdf", 
       width=5, height=4, dpi=300)


## UMAP plot by sequencing
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'sequencing',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_sequencing.pdf", 
       width=5, height=4, dpi=300)

## UMAP plot by sex
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'sex',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_sex.pdf", 
       width=5, height=4, dpi=300)


## UMAP plot by fusion
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'fusion',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_fusion.pdf", 
       width=5, height=4, dpi=300)

## UMAP plot by status
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'status',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_status.pdf", 
       width=5, height=4, dpi=300)

## UMAP plot by location
DimPlot(PDX.combined, reduction = "UMAP_no_regression", group.by = 'location',
        pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()
ggsave("/mnt/Sara/output/metadata/Patel_Danielli_Langenau/no_reg/3_UMAP_location.pdf", 
       width=5, height=4, dpi=300)


