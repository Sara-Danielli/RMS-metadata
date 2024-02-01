rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(writexl)
library(data.table)
library(viridis)
library(future)
library(paletteer)
library(SeuratObject)
library(SCpubr)
library(readxl)
library(clustree)
library(cowplot)

base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
list_dir <- '/Volumes/Sara_PhD/scRNAseq_data/list_final'

resource_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources')
source(file.path(base_dir, "codes/MANUSCRIPT_INTEGRATION/metadata/FINAL/Functions.R"))
source(file.path(resource_dir, "Plot_style_v2.R"))

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/RPCA_name/analysis')
if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}

genelist_dir <- file.path(base_dir, 'list_final')

## read file
#PDX.integrated <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_20230710.rds"))
PDX.integrated <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds"))

DefaultAssay(PDX.integrated) <- "integrated"

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
col_fusion <- c('#D3A2C2FF', '#D3A2C2FF', '#95CECFFF', 'darkblue')

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



## Export number of cells
md <- PDX.integrated@meta.data %>% as.data.table

variables <- c('cluster_names', 'seurat_clusters', 'cluster_names_aggregate',  'name', 'model', 'origin', 'sex', "sequencing", 'age', 'site', 'status', 'location','treatment')
names(variables) <- c('cluster_names', 'seurat_clusters', 'cluster_names_aggregate', 'name', 'model', 'origin', 'sex', "sequencing", 'age', 'site', 'status', 'location','treatment')

for (i in 1:length(variables)) {
  table <- md[, .N, by = c(variables[[i]])]
  write.csv(table, file.path(analysis_dir, paste0("0_Number_cells_by_", names(variables)[i], ".csv")))
}


## PC plots
variables2 <- c('name', 'model', 'origin', 'subtype')
names(variables2) <- c('name', 'model', 'origin', 'subtype')
colors_variables2 <- list(col_name, col_model, col_origin, col_subtype)

for (a in 1:length(variables2)) {
  DimPlot(PDX.integrated, reduction = "pca_rpca", group.by = variables2[[a]], cols = colors_variables2[[a]], pt.size = 2, 
          raster=TRUE, shuffle=TRUE, raster.dpi = c(1028, 1028)) + NoLegend()
  ggsave(file.path(analysis_dir, paste0("1_PCA_",names(variables2)[a],"_no_legend.pdf")), width=4, height=4, dpi=300)
}

for (a in 1:length(variables2)) {
  DimPlot(PDX.integrated, reduction = "pca_rpca", group.by = variables2[[a]], cols = colors_variables2[[a]], pt.size = 2, 
          raster=TRUE, shuffle=TRUE, raster.dpi = c(1028, 1028)) 
  ggsave(file.path(analysis_dir, paste0("1_PCA_",names(variables2)[a],".pdf")), width=10, height=4, dpi=300)
}

# PC plots of scores
scores <- c('Wei_mesenchymal1', "Wei_proliferative1", "Wei_muscle1",
            "Danielli_MuSC1",  "Danielli_cycling1", "Danielli_differentiated1",
            "Patel_mesoderm1", "Patel_myoblast1", "Patel_myocyte1",
            "Common_Stemcell1", "Common_proliferative1", "Common_differentiated1",
            "All_Stemcell1", "All_proliferative_DW1", "All_differentiated1")

p <- FeaturePlot(PDX.integrated, features = scores, reduction = "pca_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 2, raster=TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=3)
ggsave(file.path(analysis_dir, "2_PCA_scores.pdf"), width=10, height=14, dpi=300)

# PC plot of myogenic markers
myogenesis_markers <- c('CD44','CDC20', 'MYOG', 'MYH3')
p <- FeaturePlot(PDX.integrated, features = myogenesis_markers, reduction = "pca_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 1, raster=TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=4)
ggsave(file.path(analysis_dir, "3_PCA_myogenesis.pdf"), width=9, height=2, dpi=300)


## UMAP plots

# UMAP to visualize distribution of different samples/subtypes
for (a in 1:length(variables2)) {
  DimPlot(PDX.integrated, reduction = "umap_rpca", group.by = variables2[[a]], cols = colors_variables2[[a]], pt.size = 2, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
  ggsave(file.path(analysis_dir, paste0("5_UMAP_",names(variables2)[a],"_no_legend.pdf")), width=4, height=4, dpi=300)
}

for (a in 1:length(variables2)) {
  DimPlot(PDX.integrated, reduction = "umap_rpca", group.by = variables2[[a]], cols = colors_variables2[[a]], pt.size = 2, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) 
  ggsave(file.path(analysis_dir, paste0("5_UMAP_",names(variables2)[a],".pdf")), width=12, height=5, dpi=300)
}


# UMAP plot Louvain clusters 
DimPlot(PDX.integrated, reduction = "umap_rpca",  label = FALSE, cols = col_cluster_names, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, paste0("7_UMAP_clusters.pdf")), width=5, height=5, dpi=300)

DimPlot(PDX.integrated, reduction = "umap_rpca",  label = 1, cols = col_cluster_names, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, paste0("7_UMAP_clusters_labels.pdf")), width=, height=, dpi=300)

DimPlot(PDX.integrated, reduction = "umap_rpca",  label = 1, cols = col_cluster_names, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012)) + NoAxes()
ggsave(file.path(analysis_dir, paste0("7_UMAP_clusters_legend.pdf")), width=5, height=4, dpi=300)


# UMAP plots of scores
p <- FeaturePlot(PDX.integrated, features = scores, reduction = "umap_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + paletteer::scale_colour_paletteer_c("grDevices::Geyser")
}
cowplot::plot_grid(plotlist = p, ncol=3)
ggsave(file.path(analysis_dir,"8_UMAP_scores_2.pdf"), width=12, height=14, dpi=300)


# UMAP plot of cell cycle distribution
DimPlot(PDX.integrated, reduction = "umap_rpca", pt.size = 2, raster=TRUE, group.by = 'Phase', cols = c('gray76', 'black', 'grey41'), raster.dpi = c(1012, 1012), shuffle=TRUE)  + NoAxes()
ggsave(file.path(analysis_dir, "9_UMAP_cellcycle.pdf"), width=6, height=5, dpi=300)

DimPlot(PDX.integrated, reduction = "umap_rpca", pt.size = 2, raster=TRUE, group.by = 'Phase', cols = c('gray76', 'black', 'grey41'), raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend()  + NoAxes()
ggsave(file.path(analysis_dir, "9_UMAP_cellcycle_no_legend.pdf"), width=5, height=5, dpi=300)

# UMAP plot of myogenesis markers
myogenesis <- c('CD44', 'MEOX2', 'PAX7', 'CDC20', 'CDK1', 'CCNB2', 'MYOG', 'ACTC1', 'TNNT2', 'MYL1', 'MYH8', 'MYH3')
myogenesis_short <- c('CD44','CDC20', 'MYOG', 'MYH3')

p <- FeaturePlot(PDX.integrated, features = myogenesis, reduction = "umap_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=3)
ggsave(file.path(analysis_dir, paste0("10_UMAP_myogenesis.pdf")), width=9, height=10, dpi=300)


p <- FeaturePlot(PDX.integrated, features = myogenesis_short, reduction = "umap_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 1, raster=TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=4)
ggsave(file.path(analysis_dir, paste0("11_UMAP_myogenesis_short.pdf")), width=9, height=2, dpi=300)


# UMAP plot with cluster names 
DimPlot(PDX.integrated, reduction = "umap_rpca",  label = FALSE, cols = col_cluster_names, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, "7_UMAP_clusters_name.pdf"), width=4, height=4, dpi=300)

DimPlot(PDX.integrated, reduction = "umap_rpca",  label = 1, cols = col_cluster_names, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, "7_UMAP_clusters_labels_name.pdf"), width=4, height=4, dpi=300)

DimPlot(PDX.integrated, reduction = "umap_rpca",  label = 1, cols = col_cluster_names, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()
ggsave(file.path(analysis_dir, "7_UMAP_clusters_legend_name.pdf"), width=5.5, height=4, dpi=300)

# UMAP plot with cluster names (aggregate)
DimPlot(PDX.integrated, reduction = "umap_rpca",  label = 0, cols = col_cluster_names_aggregate, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE)  + NoAxes()
ggsave(file.path(analysis_dir, "7b_UMAP_clusters_aggregates.pdf"), width=5.5, height=4, dpi=300)

DimPlot(PDX.integrated, reduction = "umap_rpca",  label = 0, cols = col_cluster_names_aggregate, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, "7b_UMAP_clusters_aggregates_no_legend.pdf"), width=4, height=4, dpi=300)

# Plot UMAP Louvain clusters for each subtype
## Subset subsets based on fusion protein
eRMS <- subset(PDX.integrated, subset = fusion == 'FN-RMS')
P3F1 <- subset(PDX.integrated, subset = fusion == 'PAX3::FOXO1')
P7F1 <- subset(PDX.integrated, subset = fusion == 'PAX7::FOXO1')

## reorder clusters
P7F1$cluster_names_aggregate <- factor(x = P7F1$cluster_names_aggregate, levels = c('Progenitor','TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated', 'Differentiated', 'Apoptosis'))
P3F1$cluster_names_aggregate <- factor(x = P3F1$cluster_names_aggregate, levels = c('Progenitor','TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated', 'Differentiated', 'Apoptosis'))
eRMS$cluster_names_aggregate <- factor(x = eRMS$cluster_names_aggregate, levels = c('Progenitor','TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated', 'Differentiated', 'Apoptosis'))


## Take same amount of cells from each object          
Idents(eRMS) <- 'fusion'
Idents(P3F1) <- 'fusion'
Idents(P7F1) <- 'fusion'
eRMS <- subset(eRMS, downsample = 15000)
P3F1 <- subset(P3F1, downsample = 15000)
P7F1 <- subset(P7F1, downsample = 15000)

Idents(eRMS) <- "cluster_names_aggregate"
Idents(P3F1) <- "cluster_names_aggregate"
Idents(P7F1) <- "cluster_names_aggregate"

## Plot UMAP for each subtype 
subtype <- c(eRMS, P3F1, P7F1)
names(subtype) <- c('FNRMS', 'PAX3FOXO1', 'PAX7FOXO1')
colors_subtype <- list(col_name, col_model, col_origin, col_subtype)

for (i in 1:length(subtype)) {
  DimPlot(subtype[[i]], reduction = "umap_rpca", label = FALSE, cols = col_cluster_names_aggregate, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend()
  ggsave(file.path(analysis_dir, paste0("7c_UMAP_clusters_aggregates",names(subtype)[i],".pdf")), width=5, height=5, dpi=300)
}


## Heatmap plot of all markers subsetting to 100/200/400 cells
# Read markers
PDX.integrated.markers <- read.csv(file.path(analysis_dir, "12_markers_SCT.csv"))

order_cluster <- c('8-Quiescent Progenitor', '1-Cycling Progenitor',  '3-Cycling', 
                   '2-DNA replication', '5-DNA replication', '10-Ground', '4-Ground', '0-Ground', 
                   '11-Ground', '6-Cycling differentiated', '9-Quiescent Differentiated', '7-Apoptosis')

PDX.integrated.markers2 <- PDX.integrated.markers %>%
  mutate(cluster =  factor(cluster, levels = order_cluster)) %>%
  arrange(cluster) 

PDX.integrated.markers2 <- PDX.integrated.markers %>% group_by(cluster)

number_cells <- c(100, 200, 400)
names(number_cells) <- c('100', '200', '400')

## Heatmap plot of all markers subsetting to 100/200/400 cells
Idents(PDX.integrated) = 'seurat_clusters'
for (a in 1:length(number_cells)) {
  DoHeatmap(subset(PDX.integrated, downsample = number_cells[[a]]), features = PDX.integrated.markers2$gene, draw.lines = TRUE, angle=30, size=4, group.colors = col_cluster_names, raster=TRUE)  +
    scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')))
  ggsave(file.path(analysis_dir, paste0("8_Heatmap_down_",names(number_cells)[a],"_cells.png")), width=8, height=10, dpi=300)
}

## Heatmap plot of top10 markers subsetting to 100/200/400 cells
top10 <- PDX.integrated.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
for (a in 1:length(number_cells)) {
  DoHeatmap(subset(PDX.integrated, downsample = number_cells[[a]]), features = top10$gene, draw.lines = TRUE, angle=30, size=4, group.colors = col_cluster_names, raster=TRUE)  +
    scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')))
  ggsave(file.path(analysis_dir, paste0("8_Heatmap_down_top10",names(number_cells)[a],"_cells.png")), width=8, height=10, dpi=300)
}

## Heatmap plot of top50 markers s
top50 <- PDX.integrated.markers2 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoHeatmap(PDX.integrated, features = top50$gene, draw.lines = TRUE, angle=30, size=4, group.colors = col_cluster_names, raster=TRUE)  +
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')))
ggsave(file.path(analysis_dir, "8_Heatmap_down_top50.png"), width=8, height=10, dpi=300)


## Correlation matrix clusters
SCpubr::do_CorrelationPlot(sample = PDX.integrated, group.by = 'cluster_names', cell_size = 10)
p

## Bar plot clusters distribution
variables <- list(PDX.integrated$name, PDX.integrated$subtype, PDX.integrated$model, PDX.integrated$origin, PDX.integrated$fusion)
names(variables) <- c('name', 'subtype', 'model', 'origin', 'fusion')

# Bar plot clusters
for (a in 1:length(variables)) {
  plot_bar(PDX.integrated, variables[[a]], PDX.integrated$cluster_names, col_cluster_names) 
  ggsave(file.path(analysis_dir, paste0("9_Barplot_cluster_distribution_split_by",names(variables)[a],".pdf")), width=6, height=5, dpi=300)
}

# Bar plot cluster aggregates
for (a in 1:length(variables)) {
  plot_bar(PDX.integrated, variables[[a]], PDX.integrated$cluster_names_aggregate, col_cluster_names_aggregate) 
  ggsave(file.path(analysis_dir, paste0("9_Barplot_aggregate_cluster_distribution_split_by",names(variables)[a],".pdf")), width=6, height=5, dpi=300)
}

# Bar plot cycling property by sample
plot_bar(PDX.integrated, PDX.integrated$name, PDX.integrated$cluster_names_aggregate, col_cluster_names_aggregate) 
ggsave(file.path(analysis_dir, paste0("9_Barplot_cluster_distribution_split_by_patient.pdf")), width=16, height=4, dpi=300)

# Bar plot cycling property by cluster 
plot_bar(PDX.integrated, PDX.integrated$cluster_names, PDX.integrated$Cycling_prop, c('black', 'grey')) 
ggsave(file.path(analysis_dir, paste0("9_Barplot_cluster_distribution_split_by",names(variables)[a],".pdf")), width=6, height=4, dpi=300)




## Violin plots scores

Vln_scores <- c("Common_Stemcell1", "Common_proliferative1", "Common_differentiated1", "DS.Difference_common", "Wei_FP1", "Wei_FN1" )
names(Vln_scores) <- c("Stemcell", "Proliferative", "Differentiated", "Muscle_lineage_score", "Wei_FP", "Wei_FN" )

for (a in 1:length(Vln_scores)) {
  VlnPlot(PDX.integrated,features = Vln_scores[[a]], group.by = 'name',  split.by = 'fusion',  pt.size=0,  cols = col_fusion) 
  ggsave(file.path(analysis_dir, paste0("15_Vln_plot_scores_",names(Vln_scores)[a],".pdf")), width=15, height=4, dpi=300)
}

p1 <- VlnPlot(PDX.integrated, features = 'Wei_FP1', group.by = 'name',  split.by = 'fusion',  pt.size=0,  cols = col_fusion) 
p2 <- VlnPlot(PDX.integrated, features = 'Wei_FN1', group.by = 'name',  split.by = 'fusion',  pt.size=0,  cols = col_fusion) 

compare_means(Wei_FP1 ~ fusion, data = PDX.integrated[[]])
p3 <- VlnPlotScoresModel(P7F1.integrated, features = 'DS.Difference_common', y = 1.4) +
  stat_compare_means(method = 't.test', size = 3)

plot_grid(p1, p2, ncol = 1, align = "h", axis = "tb")
ggsave(file.path(analysis_dir,"15_VlnPlot_FP_FN_scores.pdf"), width=16, height=9, dpi=300)


p1 <- VlnPlot(PDX.integrated, features = 'Common_proliferative1', group.by = 'subtype', pt.size=0,  cols = col_subtype) 
p2 <- VlnPlot(PDX.integrated, features = 'Wei_FN1', group.by = 'name',  split.by = 'fusion',  pt.size=0,  cols = col_fusion) 



## FeaturePlot
p1 <- FeatureScatter(PDX.integrated,
               feature1 = 'DS.Difference_common',
               feature2 = "Common_proliferative1",
               group.by='name',
               shuffle=TRUE,
               cols = col_name,
               raster=TRUE,
               pt.size = 3, raster.dpi=c(1021,1012)) + NoLegend() #+ xlim(c(-1.1, 1.5)) + ylim(c(-0.4, 1.7))
ggsave(file.path(analysis_dir,"18_FeaturePlot_cell_state_plot_subtype.pdf"), width=4, height=4, dpi=300)


## FeaturePlot average
Idents(object = PDX.integrated) <- "name"
PDX.integrated2 <- DietSeurat(PDX.integrated, 
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = 'integrated',
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)
rm(PDX.integrated)

PDX.integrated <- PDX.integrated2
Idents(object = PDX.integrated) <- "name"
PDX.integrated.name.averages <- AverageExpression(PDX.integrated, assays = "integrated", return.seurat = TRUE, verbose = TRUE)

## Order samples
## Change cluster names
levels(PDX.integrated.name.averages)
PDX.integrated.name.averages <- StashIdent(PDX.integrated.name.averages, save.name = "name")

# Score each sample for cell states
Common_differentiated <- read_excel(file.path(list_dir, "Common_differentiated.xlsx"), col_names=FALSE)
Common_differentiated <- as.list(Common_differentiated)
Common_proliferative <- read_excel(file.path(list_dir, "Common_proliferative.xlsx"), col_names=FALSE)
Common_proliferative <- as.list(Common_proliferative)
Common_Stemcell <- read_excel(file.path(list_dir, "Common_Stemcell.xlsx"), col_names=FALSE)
Common_Stemcell <- as.list(Common_Stemcell)

PDX.integrated.name.averages <- AddModuleScore(object = PDX.integrated.name.averages, features = Common_differentiated, ctrl = 50, name = "Common_differentiated")
PDX.integrated.name.averages <- AddModuleScore(object = PDX.integrated.name.averages, features = Common_Stemcell, ctrl = 50, name = "Common_Stemcell")
PDX.integrated.name.averages <- AddModuleScore(object = PDX.integrated.name.averages, features = Common_proliferative, ctrl = 50, name = "Common_proliferative")
PDX.integrated.name.averages$DS.Difference_common <- PDX.integrated.name.averages$Common_differentiated1 - PDX.integrated.name.averages$Common_Stemcell1

# Scale data
ScaleData(PDX.integrated.name.averages)

# Plots
FeatureScatter(PDX.integrated.name.averages,
               feature1 = 'DS.Difference_common',
               feature2 = "Common_proliferative1",
               group.by='name',
               cols = col_name,
               pt.size = 3) 
ggsave(file.path(analysis_dir,"19_FeaturePlot_cell_state_plot_aggregate.pdf"), width=12, height=5, dpi=300)

p2 <- FeatureScatter(PDX.integrated.name.averages,
               feature1 = 'DS.Difference_common',
               feature2 = "Common_proliferative1",
               group.by='name',
               cols = col_name,
               pt.size = 3) +NoLegend()
ggsave(file.path(analysis_dir,"19_FeaturePlot_cell_state_plot_aggregate_2.pdf"), width=4, height=4, dpi=300)


p1 + p2
ggsave(file.path(analysis_dir,"20_FeaturePlot_cell_state_plot_aggregate_2.pdf"), width=8.5, height=4, dpi=300)

## Cellular states - 3 variable plot
# Upload excel list with common genes
gene_markers <- read_excel(file.path(base_dir, "gene_list.xlsx"))

gene_set <- list("Differentiated" = gene_markers$differentiated,
                 "Proliferative" = gene_markers$proliferative[1:132],
                 "Progenitor" = gene_markers$progenitor[1:171])

Idents(PDX.integrated) = 'cluster_names_aggregate'
cellular_plot_cluster2 <- SCpubr::do_CellularStatesPlot(PDX.integrated, 
                                                        input_gene_list = gene_set,
                                                        x1 = "Differentiated",
                                                        y1 = "Proliferative",
                                                        x2 = "Progenitor",
                                                        pt.size = 2,
                                                        nbin = 1,
                                                        ctrl = 100)
cellular_plot_cluster2 

ggsave(file.path(analysis_dir,"16_cell_state_plot_subtype.pdf"), width=6, height=7, dpi=300)

Idents(PDX.integrated) = 'subtype'
cellular_plot_cluster2 <- SCpubr::do_CellularStatesPlot(PDX.integrated, 
                                                        input_gene_list = gene_set,
                                                        x1 = "Differentiated",
                                                        y1 = "Proliferative",
                                                        x2 = "Progenitor",
                                                        pt.size = 2,
                                                        nbin = 1,
                                                        ctrl = 100)
cellular_plot_cluster2
ggsave(file.path(analysis_dir,"19_FeaturePlot_cell_state_plot_aggregate.pdf"), width=12, height=5, dpi=300)

 


Idents(PDX.integrated) = 'cluster_names'
PDX.integrated_small <- subset(PDX.integrated, downsample = 500)

p <- SCpubr::do_CorrelationPlot(sample = PDX.integrated_small, cell_size = 10)
ggsave(file.path(analysis_dir,"17_correlation_subset_500_cells.pdf"), width=6, height=7, dpi=300)





RMS.integrated <- readRDS(file.path(base_dir, "write/RMS_atlas_final_20240130.rds"))

# Violin plot neuronal markers
Vln_scores <- c("TBXT", "SOX2")
names(Vln_scores) <- c("TBXT", "SOX2")

for (a in 1:length(Vln_scores)) {
  VlnPlot(RMS.integrated, features = Vln_scores[[a]], group.by = 'name',  split.by = 'subtype',  pt.size=0,  cols = col_subtype) 
  ggsave(file.path(analysis_dir, paste0("18_Vln_plot_scores_",names(Vln_scores)[a],".pdf")), width=15, height=5, dpi=300)
}


# Dotplot scores
DotPlot(RMS.integrated, 
        features = names(Vln_scores), 
        group.by = 'Cluster assignment',
        split.by = 'subtype',
        assay = 'RNA', 
        cols = c("white", "red3"),
        scale = FALSE,
        col.min = 0)  + 
  #scale_colour_distiller(palette="RdBu") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title=element_blank(),
        legend.text = element_text(size = 12),
        #legend.title = element_blank()
  ) 
ggsave(file.path(plot_dir, paste0("2_DotPlot_scores_", names(seurat_objects)[i], ".pdf")), width=5.5, height=4.5, dpi=300)


# Plot Clustree  --------------------------------------------------------------------------
# re-run Louvain clusters with different resolutions
PDX.integrated <- FindClusters(PDX.integrated, resolution = seq(0, 1, 0.1))

# Plot Clustree 
clustree(PDX.integrated)
ggsave(file.path(analysis_dir,"21_Clustree_res_01_1.pdf"), width=8, height=10, dpi=300)


# Dot/Vln plot progenitor/differentiated scores for each population ------------------------------------------
DotPlot(RMS.integrated, 
        features = c("Progenitor score", "Differentiated score"), 
        group.by = 'Cluster assignment',
        assay = 'RNA', 
        cols = c("white", "red3"),
        scale = F,
        #col.min = 0
) + 
  #scale_colour_distiller(palette="RdBu") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title=element_blank(),
        legend.text = element_text(size = 12),
        #legend.title = element_blank()
  ) 
ggsave(file.path(analysis_dir, paste0("22_DotPlot_Progenitor_differentiated_score.pdf")), width=4.5, height=5, dpi=300)



scores <- c("Progenitor score", "Differentiated score")
titles <- c("Progenitor score", "Differentiated score")

p <- VlnPlot(RMS.integrated, features = scores, group.by = 'Cluster assignment', combine = FALSE, pt.size=0, cols = col_model) 

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() +
    labs (y='Module score, AU', x='', title = titles[[i]]) + 
    scale_fill_manual(values = col_cluster_names_aggregate_integrated) + 
    geom_boxplot(outlier.shape=NA, width=0.1, fill="white") + NoLegend() +
    theme_vln
}
plot_grid(plotlist = p, ncol=2)
ggsave(file.path(analysis_dir, paste0("23_VlnPlot_Progenitor_differentiated_score.pdf")), width=8, height=5, dpi=300)


