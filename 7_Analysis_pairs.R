rm(list = ls())

library(Seurat)
library(patchwork)
library(viridis)
library(future)
library(ggplot2)
library(readxl)
library(writexl)
library(data.table)
library(SingleCellExperiment)
library(paletteer)
library(SeuratObject)
library(tidyverse)
library(cowplot)
library(dplyr)

setwd('/Volumes/Sara_PhD/scRNAseq_data')
source("./codes/MANUSCRIPT_INTEGRATION/metadata/FINAL/Functions.R")


PDX.integrated <- readRDS("./write/Danielli_Patel_Langenau_20230710.rds")

## select unintegrated object
DefaultAssay(PDX.integrated) <- "RNA"

## Define colors 
    # Color model
    col_model <- c("#009B9EFF","#A7D3D4FF",  "#E4C1D9FF","#C75DABFF")
    
    # Color origin
    col_origin <- paletteer::paletteer_d("ggthemes::excel_Slice", n=4)
    
    # Color name
    col_aRMS <- paletteer::paletteer_c("ggthemes::Blue-Teal", n = 27)
    col_eRMS <- paletteer::paletteer_c("ggthemes::Purple", n = 47)
    col_name <- c(col_aRMS, col_eRMS)
    
    # Color subtype
    col_subtype <- c('#D3A2C2FF', '#95CECFFF')
    
    # Color Louvain clusters
    col_cluster_names =  paletteer::paletteer_d("rcartocolor::Pastel", n = 12)
    col_cluster_names_aggregate <- c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF", "#8E0152FF", '#D8D8E0FF')

## Order factors
    PDX.integrated$origin <- factor(x = PDX.integrated$origin, levels = c("Wei et al.", "Patel et al.", "Danielli et al.", "Weng et al."))
    PDX.integrated$model <- factor(x = PDX.integrated$model, levels = c("Patient", "O-PDX", "Primary culture", "Cell line"))
    PDX.integrated$subtype <- factor(x = PDX.integrated$subtype, levels = c("eRMS", "aRMS"))
    PDX.integrated$name <- factor(x = PDX.integrated$name, levels = c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3', 
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
    
    PDX.integrated$cluster_names <- factor(x = PDX.integrated$cluster_names, 
                                           levels = c('8-Progenitor','1-TR-progenitor', '3-Proliferative', '2-DNA replication',  '5-DNA replication', 
                                                      '10-DNA replication', '4-Ground','0-Ground', '11-Ground','6-TR-differentiated',  '9-Differentiated', '7-Apoptosis'))
    
    PDX.integrated$cluster_names_aggregate <- factor(x = PDX.integrated$cluster_names_aggregate, 
                                                     levels = c('Progenitor','TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated', 'Differentiated', 'Apoptosis'))
    
    PDX.integrated$seurat_clusters <- factor(x = PDX.integrated$seurat_clusters, levels = c('8', '1',  '3', '2', '5', '10', '4', '0', '11', '6', '9', '7'))


## Subset subtypes
    Idents(object = PDX.integrated) <- "name"
    
    mast139_pairs <- subset(PDX.integrated, idents = c( 'Mast139',  'Mast139_SC', 'SJRHB013758_X2'))
    mast111_pairs <- subset(PDX.integrated, idents = c('Mast111',  'SJRHB013758_X1'))
    mast95_pairs <- subset(PDX.integrated, idents = c('Mast95',  'SJRHB013757_X1'))
    mast39_pairs <- subset(PDX.integrated, idents = c('Mast39',  'SJRHB000026_X1'))
    mast85_pairs <- subset(PDX.integrated, idents = c('Mast85_r1', 'Mast85_r2', 'Mast85_r2_SC', 'SJRHB000026_X2'))

    # Remove integrated obejct to free up space
    rm(PDX.integrated)
    
# Perform standard downstream to susbet objects
    object_list <- list(mast139_pairs, mast111_pairs, mast95_pairs, mast39_pairs, mast85_pairs)
    names(object_list) <- c('mast139_pairs', 'mast111_pairs', 'mast95_pairs', 'mast39_pairs', 'mast85_pairs')
    
    object_list <- lapply(X = object_list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, nfeatures = 2000)
      x <- ScaleData(x,  verbose = TRUE)
      x <- RunPCA(x, verbose = TRUE)
    })

    
    # Elbow plots
    for (a in 1:length(object_list)) {
      ElbowPlot(object_list[[a]])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/0_ElbowPlot_",names(object_list)[a],".pdf"),
             width=6, height=4, dpi=300)
    }
    
    # Perform UMAP
    object_list <- lapply(X = object_list, FUN = function(x) {
      x <- FindNeighbors(x, reduction = "pca", dims = 1:15)
      x <- RunUMAP(x, dims = 1:15, reduction = "pca",  reduction.key = "umap", reduction.name = "umap")
    })
    
  
# Plot UMAP (no regression = reduction umap)
    # UMAP origin
    for (a in 1:length(object_list)) {
      DimPlot(object_list[[a]], reduction = "umap", shuffle=TRUE, group.by = 'origin', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_origin) + NoAxes() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_UMAP_",names(object_list)[a],".pdf"),width=6, height=4, dpi=300)    
    }
    
    p <- list()
    for (a in 1:length(object_list)) {
      p[[a]] <- DimPlot(object_list[[a]], reduction = "umap", shuffle=TRUE, group.by = 'origin', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_origin) + NoAxes() + NoLegend() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_UMAP_",names(object_list)[a],"_no_legend.pdf"),width=4, height=4, dpi=300)   
    }
    
    # UMAP cluster
    for (a in 1:length(object_list)) {
      DimPlot(object_list[[a]], reduction = "umap", shuffle=TRUE, group.by = 'cluster_names_aggregate', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_cluster_names_aggregate) + NoAxes() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_UMAP_cluster",names(object_list)[a],".pdf"),width=6, height=4, dpi=300)   
    }
    
    p2 <- list()
    for (a in 1:length(object_list)) {
      p2[[a]] <- DimPlot(object_list[[a]], reduction = "umap", shuffle=TRUE, group.by = 'cluster_names_aggregate', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_cluster_names_aggregate) + NoAxes()+  NoLegend() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_UMAP_cluster",names(object_list)[a],"_no_legend.pdf"),width=4, height=4, dpi=300)   
    }
    
    
    # UMAP sample name
    for (a in 1:length(object_list)) {
      DimPlot(object_list[[a]], reduction = "umap", shuffle=TRUE, group.by = 'name', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_cluster_names) + NoAxes()  + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_UMAP_sample",names(object_list)[a],".pdf"),width=6, height=4, dpi=300)   
    }
    
    p3 <- list()
    for (a in 1:length(object_list)) {
      p3[[a]] <- DimPlot(object_list[[a]], reduction = "umap", shuffle=TRUE, group.by = 'name', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_cluster_names) + NoAxes() +  NoLegend() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_UMAP_sample",names(object_list)[a],"_no_legend.pdf"),width=4, height=4, dpi=300)   
    }
    
    # UMAP model
    for (a in 1:length(object_list)) {
      DimPlot(object_list[[a]], reduction = "umap", shuffle=TRUE, group.by = 'model', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_model) + NoAxes()  + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_UMAP_model",names(object_list)[a],".pdf"),width=6, height=4, dpi=300)   
    }
    
    p4 <- list()
    for (a in 1:length(object_list)) {
      p4[[a]] <- DimPlot(object_list[[a]], reduction = "umap", shuffle=TRUE, group.by = 'model', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_model) + NoAxes() +  NoLegend() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_UMAP_model",names(object_list)[a],"_no_legend.pdf"),width=4, height=4, dpi=300)   
    }
    
    
    
    # Bar plots
    p5 <- plot_bar(mast139_pairs, mast139_pairs$origin, mast139_pairs$cluster_names_aggregate, col_cluster_names_aggregate) + NoLegend()
    p6 <- plot_bar(mast111_pairs, mast111_pairs$origin, mast111_pairs$cluster_names_aggregate, col_cluster_names_aggregate) + NoLegend()
    p7 <- plot_bar(mast95_pairs, mast95_pairs$origin, mast95_pairs$cluster_names_aggregate, col_cluster_names_aggregate) + NoLegend()
    p8 <- plot_bar(mast39_pairs, mast39_pairs$origin, mast39_pairs$cluster_names_aggregate, c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF",  '#D8D8E0FF')) + NoLegend()
    p9 <- plot_bar(mast85_pairs, mast85_pairs$origin, mast85_pairs$cluster_names_aggregate, c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF",  '#D8D8E0FF'))+ NoLegend() 
    

    
    plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], 
              p2[[1]], p2[[2]], p2[[3]], p2[[4]], p2[[5]],
              p3[[1]], p3[[2]], p3[[3]], p3[[4]], p3[[5]],
              p4[[1]], p4[[2]], p4[[3]], p4[[4]], p4[[5]],
              p5, p6, p7, p8, p9,
              ncol = 5, align = "h", axis = "tb")
    ggsave("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/3_UMAP_combined_no_legend.pdf",width=12, height=16, dpi=300)   

    
    
    
    
# Plot UMAP with RPCA correction (no regression = umap_rpca)
    # UMAP origin
    for (a in 1:length(object_list)) {
      DimPlot(object_list[[a]], reduction = "umap_rpca", shuffle=TRUE, group.by = 'origin', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_origin) + NoAxes() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_",names(object_list)[a],".pdf"),width=6, height=4, dpi=300)   
    }
    
    p <- list()
    for (a in 1:length(object_list)) {
      p[[a]] <- DimPlot(object_list[[a]], reduction = "umap_rpca", shuffle=TRUE, group.by = 'origin', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_origin) + NoAxes() + NoLegend() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_",names(object_list)[a],"_no_legend.pdf"),width=4, height=4, dpi=300)   
    }
    
    # umap_rpca cluster
    for (a in 1:length(object_list)) {
      DimPlot(object_list[[a]], reduction = "umap_rpca", shuffle=TRUE, group.by = 'cluster_names_aggregate', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_cluster_names_aggregate) + NoAxes() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_cluster",names(object_list)[a],".pdf"),width=6, height=4, dpi=300)   
    }
    
    p2 <- list()
    for (a in 1:length(object_list)) {
      p2[[a]] <- DimPlot(object_list[[a]], reduction = "umap_rpca", shuffle=TRUE, group.by = 'cluster_names_aggregate', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_cluster_names_aggregate) + NoAxes()+  NoLegend() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_cluster",names(object_list)[a],"_no_legend.pdf"),width=4, height=4, dpi=300)   
    }
    
    DimPlot(mast39_pairs, reduction = "umap_rpca", shuffle=TRUE, group.by = 'cluster_names_aggregate', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), 
            cols = c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF",  '#D8D8E0FF')) + NoAxes() + ggtitle('mast39_pairs')
    ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_clustermast39_pairs.pdf"),width=6, height=4, dpi=300)   
    
    DimPlot(mast85_pairs, reduction = "umap_rpca", shuffle=TRUE, group.by = 'cluster_names_aggregate', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), 
            cols = c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF",  '#D8D8E0FF')) + NoAxes() + ggtitle('mast85_pairs')
    ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_clustermast85_pairs.pdf"),width=6, height=4, dpi=300)   
    
    pmast39 <- DimPlot(mast39_pairs, reduction = "umap_rpca", shuffle=TRUE, group.by = 'cluster_names_aggregate', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), 
                       cols = c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF",  '#D8D8E0FF')) + NoLegend() + NoAxes() + ggtitle('mast39_pairs')
    ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_clustermast39_pairs_no_legend.pdf"),width=4, height=4, dpi=300)   
    
    pmast85 <- DimPlot(mast85_pairs, reduction = "umap_rpca", shuffle=TRUE, group.by = 'cluster_names_aggregate', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), 
                       cols = c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF",  '#D8D8E0FF')) + NoLegend()+ NoAxes()  + ggtitle('mast85_pairs')
    ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_clustermast85_pairs_no_legend.pdf"),width=4, height=4, dpi=300)   
    
    
    # umap_rpca sample name
    for (a in 1:length(object_list)) {
      DimPlot(object_list[[a]], reduction = "umap_rpca", shuffle=TRUE, group.by = 'name', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_cluster_names) + NoAxes() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_sample",names(object_list)[a],".pdf"),width=6, height=4, dpi=300)   
    }
    
    p3 <- list()
    for (a in 1:length(object_list)) {
      p3[[a]] <- DimPlot(object_list[[a]], reduction = "umap_rpca", shuffle=TRUE, group.by = 'name', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_cluster_names) + NoAxes() +  NoLegend() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_sample",names(object_list)[a],"_no_legend.pdf"),width=4, height=4, dpi=300)   
    }
    
    # umap_rpca model
    for (a in 1:length(object_list)) {
      DimPlot(object_list[[a]], reduction = "umap_rpca", shuffle=TRUE, group.by = 'model', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_model) + NoAxes()  + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_model",names(object_list)[a],".pdf"),width=6, height=4, dpi=300)   
    }
    
    p4 <- list()
    for (a in 1:length(object_list)) {
      p4[[a]] <- DimPlot(object_list[[a]], reduction = "umap_rpca", shuffle=TRUE, group.by = 'model', pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), cols = col_model) + NoAxes() +  NoLegend() + ggtitle(names(object_list)[a])
      ggsave(paste0("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/1_umap_rpca_model",names(object_list)[a],"_no_legend.pdf"),width=4, height=4, dpi=300)   
    }
    
    
    
    # Bar plots
    p5 <- plot_bar(mast139_pairs, mast139_pairs$name, mast139_pairs$cluster_names_aggregate, col_cluster_names_aggregate) + NoLegend()
    p6 <- plot_bar(mast111_pairs, mast111_pairs$name, mast111_pairs$cluster_names_aggregate, col_cluster_names_aggregate) + NoLegend()
    p7 <- plot_bar(mast95_pairs, mast95_pairs$name, mast95_pairs$cluster_names_aggregate, col_cluster_names_aggregate) + NoLegend()
    p8 <- plot_bar(mast39_pairs, mast39_pairs$name, mast39_pairs$cluster_names_aggregate, c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF",  '#D8D8E0FF')) + NoLegend()
    p9 <- plot_bar(mast85_pairs, mast85_pairs$name, mast85_pairs$cluster_names_aggregate, c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF",  '#D8D8E0FF'))+ NoLegend() 

    
    
    plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], 
              p2[[1]], p2[[2]], p2[[3]], pmast39, pmast85,
              p3[[1]], p3[[2]], p3[[3]], p3[[4]], p3[[5]],
              p4[[1]], p4[[2]], p4[[3]], p4[[4]], p4[[5]],
              p5, p6, p7, p8, p9,
              ncol = 5, align = "h", axis = "tb")
    ggsave("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/3_umap_rpca_combined_no_legend.pdf",width=12, height=16, dpi=300)   
    
    plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], 
              p2[[1]], p2[[2]], p2[[3]], pmast39, pmast85,
              p5, p6, p7, p8, p9,
              ncol = 5, align = "h", axis = "tb")
    ggsave("./output/metadata/Patel_Danielli_Langenau/RPCA_name/Pairs/4_plot_paper.pdf",width=12, height=11, dpi=300)   
    
    
    
    
    
   
    