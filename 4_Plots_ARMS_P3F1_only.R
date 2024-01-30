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
library(ggpubr)
library(cowplot)

base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
resource_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources')
source(file.path(base_dir, "codes/MANUSCRIPT_INTEGRATION/metadata/FINAL/Functions.R"))
source(file.path(resource_dir, "Plot_style_v2.R"))

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/RPCA_name/P3F1analysis')
if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}

## read file
P3F1.integrated <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_20230713.rds"))

DefaultAssay(P3F1.integrated) <- "integrated"

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
col_cluster_names =  paletteer::paletteer_d("rcartocolor::Pastel", n = 7)
col_cluster_names_aggregate <- c("#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#8E0152FF",'#66C5CCFF' , '#D8D8E0FF')

## Order factors
P3F1.integrated$origin <- factor(x = P3F1.integrated$origin, levels = c("Wei et al.", "Patel et al.", "Danielli et al.", "Weng et al."))
P3F1.integrated$model <- factor(x = P3F1.integrated$model, levels = c("Patient", "O-PDX", "Primary culture", "Cell line"))
P3F1.integrated$subtype <- factor(x = P3F1.integrated$subtype, levels = c("eRMS", "aRMS"))
P3F1.integrated$name <- factor(x = P3F1.integrated$name, levels = c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3', 
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

P3F1.integrated$cluster_names <- factor(x = P3F1.integrated$cluster_names, 
                                       levels = c('Progenitor', 'Proliferative','DNA replication', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis'))

P3F1.integrated$cluster_names_aggregate <- factor(x = P3F1.integrated$cluster_names_aggregate, 
                                                 levels = c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis'))



## PC plots
variables2 <- c('name', 'model', 'origin')
names(variables2) <- c('name', 'model', 'origin')
colors_variables2 <- list(col_name, col_model, col_origin)

for (a in 1:length(variables2)) {
  DimPlot(P3F1.integrated, reduction = "pca_rpca", group.by = variables2[[a]], cols = colors_variables2[[a]], pt.size = 4, 
          raster=TRUE, shuffle=TRUE, raster.dpi = c(1028, 1028)) + NoLegend()
  ggsave(file.path(analysis_dir, paste0("1_PCA_",names(variables2)[a],"_no_legend.pdf")), width=4, height=4, dpi=300)
}

for (a in 1:length(variables2)) {
  DimPlot(P3F1.integrated, reduction = "pca_rpca", group.by = variables2[[a]], cols = colors_variables2[[a]], pt.size = 4, 
          raster=TRUE, shuffle=TRUE, raster.dpi = c(1028, 1028)) 
  ggsave(file.path(analysis_dir, paste0("1_PCA_",names(variables2)[a],".pdf")), width=10, height=4, dpi=300)
}

# PC plots of scores
scores <- c('Wei_mesenchymal1', "Wei_proliferative1", "Wei_muscle1",
            "Danielli_MuSC1",  "Danielli_cycling1", "Danielli_differentiated1",
            "Patel_mesoderm1", "Patel_myoblast1", "Patel_myocyte1",
            "Common_Stemcell1", "Common_proliferative1", "Common_differentiated1",
            "All_Stemcell1", "All_proliferative_DW1", "All_differentiated1")

p <- FeaturePlot(P3F1.integrated, features = scores, reduction = "pca_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 4, raster=TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=3)
ggsave(file.path(analysis_dir, "2_PCA_scores.pdf"), width=10, height=3, dpi=300)

# PC plot of myogenic markers
myogenesis_markers <- c('CD44','CDC20', 'MYOG', 'MYH3')
p <- FeaturePlot(P3F1.integrated, features = myogenesis_markers, reduction = "pca_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 4, raster=TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=4)
ggsave(file.path(analysis_dir, "3_PCA_myogenesis.pdf"), width=9, height=2, dpi=300)


## UMAP plots

# UMAP to visualize distribution of different samples/subtypes
for (a in 1:length(variables2)) {
  DimPlot(P3F1.integrated, reduction = "umap_rpca", group.by = variables2[[a]], cols = colors_variables2[[a]], pt.size = 4, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
  ggsave(file.path(analysis_dir, paste0("5_UMAP_",names(variables2)[a],"_no_legend.pdf")), width=4, height=4, dpi=300)
}

for (a in 1:length(variables2)) {
  DimPlot(P3F1.integrated, reduction = "umap_rpca", group.by = variables2[[a]], cols = colors_variables2[[a]], pt.size = 4, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) 
  ggsave(file.path(analysis_dir, paste0("5_UMAP_",names(variables2)[a],".pdf")), width=12, height=5, dpi=300)
}


# UMAP plot Louvain clusters 
Idents(P3F1.integrated) = P3F1.integrated$cluster_names
DimPlot(P3F1.integrated, reduction = "umap_rpca",  label = FALSE, cols = col_cluster_names, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, paste0("7_UMAP_clusters.pdf")), width=5, height=5, dpi=300)

DimPlot(P3F1.integrated, reduction = "umap_rpca",  label = 1, cols = col_cluster_names, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, paste0("7_UMAP_clusters_labels.pdf")), width=, height=, dpi=300)

DimPlot(P3F1.integrated, reduction = "umap_rpca",  label = 1, cols = col_cluster_names, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012)) + NoAxes()
ggsave(file.path(analysis_dir, paste0("7_UMAP_clusters_legend.pdf")), width=5, height=4, dpi=300)


# UMAP plots of scores
p <- FeaturePlot(P3F1.integrated, features = scores, reduction = "umap_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + paletteer::scale_colour_paletteer_c("grDevices::Geyser")
}
cowplot::plot_grid(plotlist = p, ncol=3)
ggsave(file.path(analysis_dir,"8_UMAP_scores_2.pdf"), width=10, height=3, dpi=300)


# UMAP plot of cell cycle distribution
DimPlot(P3F1.integrated, reduction = "umap_rpca", pt.size = 4, raster=TRUE, group.by = 'Phase', cols = c('black', 'gray76', 'grey41'), raster.dpi = c(1012, 1012), shuffle=TRUE)  + NoAxes()
ggsave(file.path(analysis_dir, "9_UMAP_cellcycle.pdf"), width=6, height=5, dpi=300)

DimPlot(P3F1.integrated, reduction = "umap_rpca", pt.size = 4, raster=TRUE, group.by = 'Phase', cols = c('black', 'gray76', 'grey41'), raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend()  + NoAxes()
ggsave(file.path(analysis_dir, "9_UMAP_cellcycle_no_legend.pdf"), width=5, height=5, dpi=300)

# UMAP plot of myogenesis markers
myogenesis <- c('CD44', 'MEOX2', 'PAX7', 'CDC20', 'CDK1', 'CCNB2', 'MYOG', 'ACTC1', 'TNNT2', 'MYL1', 'MYH8', 'MYH3')
myogenesis_short <- c('CD44','CDC20', 'MYOG', 'MYH3')

p <- FeaturePlot(P3F1.integrated, features = myogenesis, reduction = "umap_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=3)
ggsave(file.path(analysis_dir, paste0("10_UMAP_myogenesis.pdf")), width=9, height=10, dpi=300)


p <- FeaturePlot(P3F1.integrated, features = myogenesis_short, reduction = "umap_rpca", combine = FALSE, sort.cell = TRUE, pt.size = 1, raster=TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=4)
ggsave(file.path(analysis_dir, paste0("11_UMAP_myogenesis_short.pdf")), width=9, height=2, dpi=300)


# UMAP plot with cluster names 
Idents(P3F1.integrated) = 'cluster_names'
DimPlot(P3F1.integrated, reduction = "umap_rpca",  label = FALSE, cols = col_cluster_names, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, "7_UMAP_clusters_name.pdf"), width=4, height=4, dpi=300)

DimPlot(P3F1.integrated, reduction = "umap_rpca",  label = 1, cols = col_cluster_names, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, "7_UMAP_clusters_labels_name.pdf"), width=4, height=4, dpi=300)

DimPlot(P3F1.integrated, reduction = "umap_rpca",  label = 1, cols = col_cluster_names, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoAxes()
ggsave(file.path(analysis_dir, "7_UMAP_clusters_legend_name.pdf"), width=5.5, height=4, dpi=300)

# UMAP plot with cluster names (aggregate)
Idents(P3F1.integrated) = 'cluster_names_aggregate'
DimPlot(P3F1.integrated, reduction = "umap_rpca",  label = 0, cols = col_cluster_names_aggregate, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE)  + NoAxes()
ggsave(file.path(analysis_dir, "7b_UMAP_clusters_aggregates.pdf"), width=5.5, height=4, dpi=300)

DimPlot(P3F1.integrated, reduction = "umap_rpca",label = 0, cols = col_cluster_names_aggregate, pt.size = 4, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, "7b_UMAP_clusters_aggregates_no_legend.pdf"), width=4, height=4, dpi=300)

P3F1.integrated$model <- factor(x = P3F1.integrated$model, levels = c("Cell line", "Primary culture", "O-PDX", "Patient"))
DimPlot(P3F1.integrated, reduction = "umap_rpca", split.by = 'model',
        label = 0, cols = col_cluster_names_aggregate, pt.size = 6, raster=TRUE, raster.dpi = c(1012, 1012), shuffle=TRUE) + NoLegend() + NoAxes()
ggsave(file.path(analysis_dir, "7c_UMAP_clusters_aggregates_split_model.pdf"), width=12, height=3, dpi=300)

P3F1.integrated$model <- factor(x = P3F1.integrated$model, levels = c("Patient", "O-PDX", "Primary culture", "Cell line"))


## Heatmap plot of all markers subsetting to 100/200/400 cells
# Read markers
P3F1.integrated.markers <- read.csv(file.path(analysis_dir, "12_markers_SCT.csv"))

order_cluster <- c('Progenitor', 'Proliferative','DNA replication', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')

P3F1.integrated.markers2 <- P3F1.integrated.markers %>%
  mutate(cluster =  factor(cluster, levels = order_cluster)) %>%
  arrange(cluster) 

P3F1.integrated.markers2 <- P3F1.integrated.markers %>% group_by(cluster)

number_cells <- c(100, 200, 400)
names(number_cells) <- c('100', '200', '400')

## Heatmap plot of all markers subsetting to 100/200/400 cells
Idents(P3F1.integrated) = 'seurat_clusters'
for (a in 1:length(number_cells)) {
  DoHeatmap(subset(P3F1.integrated, downsample = number_cells[[a]]), features = P3F1.integrated.markers2$gene, draw.lines = TRUE, angle=30, size=4, group.colors = col_cluster_names, raster=TRUE)  +
    scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')))
  ggsave(file.path(analysis_dir, paste0("8_Heatmap_down_",names(number_cells)[a],"_cells.png")), width=8, height=10, dpi=300)
}

## Heatmap plot of top10 markers subsetting to 100/200/400 cells
top10 <- P3F1.integrated.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
for (a in 1:length(number_cells)) {
  DoHeatmap(subset(P3F1.integrated, downsample = number_cells[[a]]), features = top10$gene, draw.lines = TRUE, angle=30, size=4, group.colors = col_cluster_names, raster=TRUE)  +
    scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')))
  ggsave(file.path(analysis_dir, paste0("8_Heatmap_down_top10",names(number_cells)[a],"_cells.png")), width=8, height=10, dpi=300)
}

## Heatmap plot of top50 markers s
top50 <- P3F1.integrated.markers2 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoHeatmap(P3F1.integrated, features = top50$gene, draw.lines = TRUE, angle=30, size=4, group.colors = col_cluster_names, raster=TRUE)  +
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')))
ggsave(file.path(analysis_dir, "8_Heatmap_down_top50.png"), width=8, height=10, dpi=300)


## Correlation matrix clusters
SCpubr::do_CorrelationPlot(sample = P3F1.integrated, group.by = 'cluster_names', cell_size = 10)
p

## Bar plot clusters distribution
variables <- list(P3F1.integrated$name, P3F1.integrated$subtype, P3F1.integrated$model, P3F1.integrated$origin, P3F1.integrated$fusion)
names(variables) <- c('name', 'model', 'origin', 'fusion')

# Bar plot clusters
for (a in 1:length(variables)) {
  plot_bar(P3F1.integrated, variables[[a]], P3F1.integrated$cluster_names, col_cluster_names) 
  ggsave(file.path(analysis_dir, paste0("9_Barplot_cluster_distribution_split_by",names(variables)[a],".pdf")), width=6, height=5, dpi=300)
}

# Bar plot cluster aggregates
for (a in 1:length(variables)) {
  plot_bar(P3F1.integrated, variables[[a]], P3F1.integrated$cluster_names_aggregate, col_cluster_names_aggregate) 
  ggsave(file.path(analysis_dir, paste0("9_Barplot_aggregate_cluster_distribution_split_by",names(variables)[a],".pdf")), width=6, height=4, dpi=300)
}

# Bar plot cycling property by sample
plot_bar(P3F1.integrated, P3F1.integrated$name, P3F1.integrated$cluster_names_aggregate, col_cluster_names_aggregate) 
ggsave(file.path(analysis_dir, paste0("9_Barplot_cluster_distribution_split_by_patient.pdf")), width=16, height=4, dpi=300)

# Bar plot cycling property by cluster 
plot_bar(P3F1.integrated, P3F1.integrated$cluster_names, P3F1.integrated$Cycling_prop, c('black', 'grey')) 
ggsave(file.path(analysis_dir, paste0("9_Barplot_cluster_distribution_split_by",names(variables)[a],".pdf")), width=6, height=4, dpi=300)




## Violin plots scores
Vln_scores <- c("Common_Stemcell1", "Common_proliferative1", "Common_differentiated1", "DS.Difference_common")
names(Vln_scores) <- c("Stemcell", "Proliferative", "Differentiated", "Muscle_lineage_scoe")

for (a in 1:length(Vln_scores)) {
  VlnPlot(P3F1.integrated,features = Vln_scores[[a]], group.by = 'name',  split.by = 'subtype',  pt.size=0,  cols = col_subtype) 
  ggsave(file.path(analysis_dir, paste0("15_Vln_plot_scores_",names(Vln_scores)[a],".pdf")), width=15, height=5, dpi=300)
}


## FeaturePlot
FeatureScatter(P3F1.integrated,
               feature1 = 'DS.Difference_common',
               feature2 = "Common_proliferative1",
               group.by='name',
               shuffle=TRUE,
               cols = col_name,
               raster=TRUE,
               pt.size = 3, raster.dpi=c(1021,1012)) + NoLegend() #+ xlim(c(-1.1, 1.5)) + ylim(c(-0.4, 1.7))
ggsave(file.path(analysis_dir,"18_FeaturePlot_cell_state_plot_subtype.pdf"), width=6, height=7, dpi=300)


Idents(P3F1.integrated) = 'cluster_names'
P3F1.integrated_small <- subset(P3F1.integrated, downsample = 500)

SCpubr::do_CorrelationPlot(sample = P3F1.integrated_small, cell_size = 10)


## Violin plot of muscle lineage and proliferation scores across different models

P3F1.integrated_number_cells <- P3F1.integrated@meta.data %>% 
  group_by(P3F1.integrated@meta.data$model, P3F1.integrated@meta.data$name) %>% 
  summarise(n = n()) 


compare_means(DS.Difference_common ~ model, data = P3F1.integrated[[]])
my_comparisons <- list(c('Patient', 'O-PDX'), c('Patient', 'Primary culture'), c('Patient', 'Cell line'), c('Patient', 'O-PDX'),
                       c('O-PDX', 'Primary culture'), c('Primary culture', 'Cell line'), c('O-PDX', 'Cell line'))

p1 <- VlnPlotScoresModel(P3F1.integrated, features = 'DS.Difference_common', y = 3) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.adj', size = 3) 
p2 <- VlnPlotScoresModel(P3F1.integrated, features = 'Common_proliferative1', y = 3) +
  stat_compare_means(comparisons = my_comparisons, label = 'p.adj', size = 3) 

plot_grid(p2, p1, ncol = 2, align = "h", axis = "tb")
ggsave(file.path(analysis_dir,"19_VlnPlot_statistics.pdf"), width=8, height=7, dpi=300)


p1 <- VlnPlotScoresModel(P3F1.integrated, features = 'DS.Difference_common', y = 1.2) 
p2 <- VlnPlotScoresModel(P3F1.integrated, features = 'Common_proliferative1', y = 1.2) 

plot_grid(p2, p1, ncol = 2, align = "h", axis = "tb")
ggsave(file.path(analysis_dir,"19_VlnPlot.pdf"), width=6, height=4, dpi=300)




P3F1 <- readRDS(file.path(base_dir, "write/FPRMS_PAX3FOXO1_final_20240130.rds"))

# Violin plot neuronal markers
Vln_scores <- c("TBXT", "SOX2")
names(Vln_scores) <- c("TBXT", "SOX2")

for (a in 1:length(Vln_scores)) {
  VlnPlot(P3F1, features = Vln_scores[[a]], group.by = 'Cluster assignment', pt.size=0,  cols = col_cluster_names_aggregate) 
  ggsave(file.path(analysis_dir, paste0("20_Vln_plot_scores_",names(Vln_scores)[a],".pdf")), width=5, height=5, dpi=300)
}


# Dotplot scores
DotPlot(P3F1, 
        features = names(Vln_scores), 
        group.by = 'Cluster assignment',
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
ggsave(file.path(analysis_dir, paste0("21_DotPlot_scores_TBXT_SOX2.pdf")), width=5.5, height=4.5, dpi=300)
