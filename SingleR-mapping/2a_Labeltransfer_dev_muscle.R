rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratObject)
library(SingleCellExperiment)
library(readxl)
library(data.table)
library(paletteer)
library(RColorBrewer)
library(scater)
library(scRNAseq)
library(SingleR)

setwd('/Volumes/Sara_PhD/scRNAseq_data/')

#################################################
############### P3F1 FP-RMS TUMORS ###############
#################################################


######### load data ##########    
      FPRMS <- readRDS("./write/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_annotated.rds")
      muscle <- readRDS("./write/Dev_muscle_Xietal.rds")
      
      
      # reorder labels
      muscle$orig.ident.short_cell.type <- factor(x = muscle$orig.ident.short_cell.type, 
                                                         levels = c('Wk5.0_AP32_MP',  'Wk5.0_AP37_MP',  
                                                                    'Wk6.0_AP5_MP', 'Wk6.0_AP13_MP', 'Wk6.5_AP22_MB-MC', 'Wk6.5_AP22_MP',
                                                                    'Wk7.0_AP24_MB-MC','Wk7.0_AP24_MP',  'Wk7.0_AP30_MB-MC', 
                                                                    'Wk7.0_AP30_MP',  'Wk7.0_AP6_MB-MC', 'Wk7.0_AP6_MP', 
                                                                    'Wk7.25_AP29_MB',  'Wk7.25_AP29_MC', 'Wk7.25_AP29_MP', 'Wk7.25_AP29_SkM.Mesen', 'Wk7.25_AP38_MB', 
                                                                    'Wk7.25_AP38_MC',  'Wk7.25_AP38_MP','Wk7.25_AP38_SkM.Mesen',
                                                                    'Wk7.75_AP10_MB',   'Wk7.75_AP10_MC',   'Wk7.75_AP10_MP', 'Wk7.75_AP10_SkM.Mesen', 
                                                                    'Wk9_AP14_MB',  'Wk9_AP14_MC', 'Wk9_AP14_MP','Wk9_AP14_SkM.Mesen','Wk9_AP36_MB', 
                                                                    'Wk9_AP36_MC', 'Wk9_AP36_SkM.Mesen', 'Wk9_AP36_MP',
                                                                    'Wk12_AP23_MB-MC',  'Wk12_AP23_MP',  'Wk12_AP23_SkM.Mesen', 
                                                                    'Wk14_AP25_MB-MC',   'Wk14_AP25_SkM.Mesen', 'Wk14_AP25_MP',
                                                                    'Wk17_AP17_MB-MC', 'Wk17_AP17_MP', 'Wk17_AP17_SkM.Mesen',
                                                                    'Wk18_AP16_MB-MC', 'Wk18_AP16_MP', 'Wk18_AP16_SkM.Mesen',
                                                                    'Yr7_AP11_SC', 'Yr11_AP34_SC', 'Yr34_AP20_SC', 'Yr42_AP28_SC'))
      
      
# transform into single cell object
      muscle.sce <- as.SingleCellExperiment(muscle)
      FPRMS.sce <- as.SingleCellExperiment(FPRMS)

      
############################################################
      ## transfer muscle labels onto P3F1 FP-RMS tumors
############################################################
      
# use SingleR to transfer labels
      pred.FPRMS <- SingleR(test=FPRMS.sce, ref=muscle.sce, labels=muscle.sce$orig.ident.short_cell.type, de.method="wilcox")
      table <- table(pred.FPRMS$labels)
      write.csv(table, "./output/metadata/dev_muscle/SingleR_projection_FPRMS/0_prediction.csv")
      
#plot heatmap
      pdf(file="./output/metadata/dev_muscle/SingleR_projection_FPRMS/1_heatmap.pdf",
          width=9, height=5)
      plotScoreHeatmap(pred.FPRMS, order.by='label')
      dev.off()
      
# plot  per-cell ???deltas??? = difference between the score for the assigned label and the median across all labels for each cell.
      pdf(file="./output/metadata/dev_muscle/SingleR_projection_FPRMS/2_deltas.pdf",
          width=6, height=10)
      plotDeltaDistribution(pred.FPRMS, ncol = 3)
      dev.off()
      
# add predicted labels to Seurat object
FPRMS[["SingleR.labels"]] <- pred.FPRMS$labels



# reorder labels
FPRMS$SingleR.labels <- factor(x = FPRMS$SingleR.labels , 
                                                   levels = c('Wk5.0_AP32_MP',  'Wk5.0_AP37_MP',  
                                                              'Wk6.0_AP5_MP', 'Wk6.0_AP13_MP', 'Wk6.5_AP22_MB-MC', 'Wk6.5_AP22_MP',
                                                              'Wk7.0_AP24_MB-MC','Wk7.0_AP24_MP',  'Wk7.0_AP30_MB-MC', 
                                                              'Wk7.0_AP30_MP',  'Wk7.0_AP6_MB-MC', 'Wk7.0_AP6_MP', 
                                                              'Wk7.25_AP29_MB',  'Wk7.25_AP29_MC', 'Wk7.25_AP29_MP', 'Wk7.25_AP29_SkM.Mesen', 'Wk7.25_AP38_MB', 
                                                              'Wk7.25_AP38_MC',  'Wk7.25_AP38_MP',
                                                              'Wk7.75_AP10_MB',   'Wk7.75_AP10_MC',   'Wk7.75_AP10_MP', 'Wk7.75_AP10_SkM.Mesen', 
                                                              'Wk9_AP14_MB',  'Wk9_AP14_MC', 'Wk9_AP14_MP','Wk9_AP14_SkM.Mesen','Wk9_AP36_MB', 
                                                              'Wk9_AP36_MC', 'Wk9_AP36_SkM.Mesen', 'Wk9_AP36_MP',
                                                              'Wk12_AP23_MB-MC',  'Wk12_AP23_MP',  'Wk12_AP23_SkM.Mesen', 
                                                              'Wk14_AP25_MB-MC',   'Wk14_AP25_SkM.Mesen', 'Wk14_AP25_MP',
                                                              'Wk17_AP17_MB-MC', 
                                                              'Wk18_AP16_MB-MC',  'Wk18_AP16_SkM.Mesen'))


#colors time points
col1 <- c('#9E0142FF', '#9E0142FF', 
          '#D53E4FFF', '#D53E4FFF', '#D53E4FFF', '#D53E4FFF',
          '#F46D43FF', '#F46D43FF','#F46D43FF',
          '#F46D43FF','#F46D43FF','#F46D43FF',
          '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  
          '#FDAE61FF',  '#FDAE61FF',  
          '#FEE08BFF','#FEE08BFF','#FEE08BFF','#FEE08BFF',
          '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  
          '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF', 
          '#E6F598FF', '#E6F598FF', '#E6F598FF', 
          '#ABDDA4FF',  '#ABDDA4FF', '#ABDDA4FF', 
          '#66C2A5FF',  
          '#3288BDFF',   '#3288BDFF')


#colors populations
col2 <- c("#C70E7BFF", "#C70E7BFF", 
          "#C70E7BFF", "#C70E7BFF", '#FC6882FF', "#C70E7BFF",
          '#FC6882FF', "#C70E7BFF", '#FC6882FF', 
          "#C70E7BFF", '#FC6882FF', "#C70E7BFF",
          '#A6E000FF', '#1BB6AFFF', "#C70E7BFF", '#6C6C9DFF', '#A6E000FF',
          '#1BB6AFFF', "#C70E7BFF",
          '#A6E000FF','#1BB6AFFF', "#C70E7BFF", '#6C6C9DFF', 
          '#A6E000FF','#1BB6AFFF',  "#C70E7BFF", '#6C6C9DFF','#A6E000FF',
          '#1BB6AFFF', '#6C6C9DFF', "#C70E7BFF", 
          '#FC6882FF', "#C70E7BFF", '#6C6C9DFF',
          '#FC6882FF', '#6C6C9DFF', '#C70E7BFF', 
          '#FC6882FF',
          '#FC6882FF','#6C6C9DFF')


DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, label=TRUE, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS/3_UMAP_projections.pdf", width=4, height=4, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, label=TRUE, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS/3_UMAP_projections_pop.pdf", width=4, height=4, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS/3_UMAP_projections_noleg.pdf", width=4, height=4, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS/3_UMAP_projections_pop_noleg.pdf", width=4, height=4, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012))
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS/3_UMAP_projections_legend.pdf", width=9, height=5, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012))
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS/3_UMAP_projections_pop_legend.pdf", width=9, height=6, dpi=300)



#Bar plot clusters
# for each sample
plot_bar <-ggplot(FPRMS@meta.data, aes(x=name, fill= FPRMS@meta.data$SingleR.labels)) + 
  geom_bar(position = "fill", color="black") +
  scale_fill_manual(values = col1) +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_bar
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS/4_cell_composition.pdf", 
       width=10, height=6, dpi=300)

plot_bar <-ggplot(FPRMS@meta.data, aes(x=name, fill= FPRMS@meta.data$SingleR.labels)) + 
  geom_bar(position = "fill", color="black") +
  scale_fill_manual(values = col2) +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_bar
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS/4_cell_composition_pop2.pdf", 
       width=10, height=6, dpi=300)

 
## Calculate number of cells within each cluster
library(data.table)
md <- FPRMS@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "./output/metadata/dev_muscle/SingleR_projection_FPRMS/5_Cluster_information_name.csv")



#save file
saveRDS(FPRMS, file = "./write/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_annotated_label.rds")
FPRMS <- readRDS("./write/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_annotated_label.rds")










#################################################
############### P7F1 FP-RMS TUMORS ###############
#################################################


######### load data ##########    
FPRMS <- readRDS("./write/Danielli_Patel_Langenau_RPCA_ARMS_P7F1_annotated.rds")
muscle <- readRDS("./write/Dev_muscle_Xietal.rds")

# transform into single cell object
muscle.sce <- as.SingleCellExperiment(muscle)
FPRMS.sce <- as.SingleCellExperiment(FPRMS)


############################################################
## transfer muscle labels onto P3F1 FP-RMS tumors
############################################################

# use SingleR to transfer labels
pred.FPRMS <- SingleR(test=FPRMS.sce, ref=muscle.sce, labels=muscle.sce$orig.ident.short_cell.type, de.method="wilcox", assay.type.test=1)
table <- table(pred.FPRMS$labels)
write.csv(table, "./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/0_prediction.csv")

#plot heatmap
pdf(file="./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/1_heatmap.pdf",
    width=9, height=5)
plotScoreHeatmap(pred.FPRMS, order.by='label')
dev.off()

# plot  per-cell ???deltas??? = difference between the score for the assigned label and the median across all labels for each cell.
pdf(file="./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/2_deltas.pdf",
    width=6, height=10)
plotDeltaDistribution(pred.FPRMS, ncol = 3)
dev.off()

# add predicted labels to Seurat object
FPRMS[["SingleR.labels"]] <- pred.FPRMS$labels

# reorder labels
FPRMS$SingleR.labels <- factor(x = FPRMS$SingleR.labels, 
                               levels = c('Wk5.0_AP32_MP',  'Wk5.0_AP37_MP',  
                                          'Wk6.0_AP5_MP',  'Wk6.0_AP13_MP', 'Wk6.5_AP22_MB-MC', 'Wk6.5_AP22_MP',
                                          'Wk7.0_AP24_MB-MC','Wk7.0_AP24_MP',  'Wk7.0_AP30_MB-MC', 'Wk7.0_AP30_MP', 'Wk7.0_AP6_MP', 
                                          'Wk7.25_AP29_MB',  'Wk7.25_AP29_MC', 'Wk7.25_AP29_MP', 'Wk7.25_AP38_MB',  'Wk7.25_AP38_MC',  'Wk7.25_AP38_MP',  
                                          'Wk7.75_AP10_MB',   'Wk7.75_AP10_MC',   'Wk7.75_AP10_MP', 'Wk7.75_AP10_SkM.Mesen', 
                                          'Wk9_AP14_MB',  'Wk9_AP14_MC', 'Wk9_AP14_MP','Wk9_AP14_SkM.Mesen','Wk9_AP36_MB', 'Wk9_AP36_MC', 'Wk9_AP36_SkM.Mesen', 'Wk9_AP36_MP',
                                          'Wk12_AP23_MB-MC',  'Wk12_AP23_MP',  'Wk12_AP23_SkM.Mesen', 
                                          'Wk14_AP25_MB-MC',   'Wk14_AP25_SkM.Mesen', 'Wk14_AP25_MP',
                                          'Wk17_AP17_MB-MC', 
                                          'Wk18_AP16_MB-MC', 'Wk18_AP16_MP', 'Wk18_AP16_SkM.Mesen', 
                                          'Yr11_AP34_SC'))

#colors populations
col2 <- c("#C70E7BFF", "#C70E7BFF", 
          "#C70E7BFF", "#C70E7BFF", '#FC6882FF', "#C70E7BFF",
          '#FC6882FF', "#C70E7BFF", '#FC6882FF', 
          "#C70E7BFF", "#C70E7BFF",
          '#A6E000FF', '#1BB6AFFF', "#C70E7BFF",  '#6C6C9DFF',
          '#1BB6AFFF', "#C70E7BFF", 
          '#A6E000FF','#1BB6AFFF', "#C70E7BFF", '#6C6C9DFF', 
          '#A6E000FF','#1BB6AFFF',  "#C70E7BFF", '#6C6C9DFF','#A6E000FF',
          '#1BB6AFFF', '#6C6C9DFF', "#C70E7BFF", 
          '#FC6882FF', "#C70E7BFF", '#6C6C9DFF',
          '#FC6882FF', '#6C6C9DFF', "#C70E7BFF", 
          '#FC6882FF', 
          '#FC6882FF', "#C70E7BFF", '#6C6C9DFF', 
          '#172869FF')



#colors time points
col1 <- c('#9E0142FF', '#9E0142FF',
          '#D53E4FFF',  '#D53E4FFF', '#D53E4FFF', '#D53E4FFF',
          '#F46D43FF', '#F46D43FF','#F46D43FF', '#F46D43FF', '#F46D43FF',
          '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF', 
          '#FEE08BFF',  '#FEE08BFF', '#FEE08BFF', '#FEE08BFF',
          '#FFFFBFFF',  '#FFFFBFFF',   '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF', 
          '#E6F598FF',     '#E6F598FF', '#E6F598FF', 
          '#ABDDA4FF',  '#ABDDA4FF',  '#ABDDA4FF', 
          '#66C2A5FF', 
          '#3288BDFF',  '#3288BDFF', '#3288BDFF',
          '#5E4FA2FF')


DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, label=TRUE, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/3_UMAP_projections.pdf", width=4, height=4, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, label=TRUE, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/3_UMAP_projections_pop.pdf", width=4, height=4, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/3_UMAP_projections_noleg.pdf", width=4, height=4, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/3_UMAP_projections_pop_noleg.pdf", width=4, height=4, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012))
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/3_UMAP_projections_legend.pdf", width=9, height=5, dpi=300)

DimPlot(FPRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012))
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/3_UMAP_projections_pop_legend.pdf", width=9, height=5, dpi=300)



#Bar plot clusters
# for each sample
plot_bar <-ggplot(FPRMS@meta.data, aes(x=name, fill= FPRMS@meta.data$SingleR.labels)) + 
  geom_bar(position = "fill", color="black") +
  scale_fill_manual(values = col1) +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_bar
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/4_cell_composition.pdf", 
       width=10, height=6, dpi=300)

plot_bar <-ggplot(FPRMS@meta.data, aes(x=name, fill= FPRMS@meta.data$SingleR.labels)) + 
  geom_bar(position = "fill", color="black") +
  scale_fill_manual(values = col2) +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_bar
ggsave("./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/4_cell_composition_pop.pdf", 
       width=10, height=6, dpi=300)


## Calculate number of cells within each cluster
library(data.table)
md <- FPRMS@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "./output/metadata/dev_muscle/SingleR_projection_FPRMS_P7F1/5_Cluster_information_name.csv")



#save file
saveRDS(FPRMS, file = "./write/Danielli_Patel_Langenau_RPCA_ARMS_P7F1_annotated_label.rds")
FPRMS <- readRDS("./write/Danielli_Patel_Langenau_RPCA_ARMS_P7F1_annotated_label.rds")












#################################################
############### FN-RMS TUMORS ###############
#################################################


######### load data ##########    
FNRMS <- readRDS("./write/Danielli_Patel_Langenau_RPCA_ERMS_annotated.rds")
muscle <- readRDS("./write/Dev_muscle_Xietal.rds")

# transform into single cell object
muscle.sce <- as.SingleCellExperiment(muscle)
FNRMS.sce <- as.SingleCellExperiment(FNRMS)


############################################################
## transfer muscle labels onto P3F1 FP-RMS tumors
############################################################

# use SingleR to transfer labels
pred.FNRMS <- SingleR(test=FNRMS.sce, ref=muscle.sce, labels=muscle.sce$orig.ident.short_cell.type, de.method="wilcox", assay.type.test=1)
table <- table(pred.FNRMS$labels)
write.csv(table, "./output/metadata/dev_muscle/SingleR_projection_FNRMS/0_prediction.csv")

#plot heatmap
pdf(file="./output/metadata/dev_muscle/SingleR_projection_FNRMS/1_heatmap.pdf",
    width=9, height=5)
plotScoreHeatmap(pred.FNRMS, order.by='label')
dev.off()

# plot  per-cell ???deltas??? = difference between the score for the assigned label and the median across all labels for each cell.
pdf(file="./output/metadata/dev_muscle/SingleR_projection_FNRMS/2_deltas.pdf",
    width=6, height=10)
plotDeltaDistribution(pred.FNRMS, ncol = 3)
dev.off()

# add predicted labels to Seurat object
FNRMS[["SingleR.labels"]] <- pred.FNRMS$labels

# reorder labels
FNRMS$SingleR.labels <- factor(x = FNRMS$SingleR.labels, 
                               levels = c('Wk5.0_AP32_MP',  'Wk5.0_AP37_MP',  
                                          'Wk6.0_AP5_MP', 'Wk6.0_AP13_MP', 'Wk6.5_AP22_MB-MC', 'Wk6.5_AP22_MP',
                                          'Wk7.0_AP24_MB-MC','Wk7.0_AP24_MP',  'Wk7.0_AP30_MB-MC', 
                                          'Wk7.0_AP30_MP',  'Wk7.0_AP6_MB-MC', 'Wk7.0_AP6_MP', 
                                          'Wk7.25_AP29_MB',  'Wk7.25_AP29_MC', 'Wk7.25_AP29_MP', 'Wk7.25_AP29_SkM.Mesen', 'Wk7.25_AP38_MB', 
                                          'Wk7.25_AP38_MC',  'Wk7.25_AP38_MP',
                                          'Wk7.75_AP10_MB',   'Wk7.75_AP10_MC',   'Wk7.75_AP10_MP', 'Wk7.75_AP10_SkM.Mesen', 
                                          'Wk9_AP14_MB',  'Wk9_AP14_MC', 'Wk9_AP14_MP','Wk9_AP14_SkM.Mesen','Wk9_AP36_MB', 
                                          'Wk9_AP36_MC', 'Wk9_AP36_SkM.Mesen', 'Wk9_AP36_MP', 
                                          'Wk12_AP23_MB-MC',  'Wk12_AP23_MP',  'Wk12_AP23_SkM.Mesen', 
                                          'Wk14_AP25_SkM.Mesen', 'Wk14_AP25_MP',
                                          'Wk17_AP17_SkM.Mesen',
                                         'Wk18_AP16_SkM.Mesen',  'Wk18_AP16_MP',
                                         'Yr42_AP28_SC'))

#colors populations
col2 <- c("#C70E7BFF", "#C70E7BFF", 
          "#C70E7BFF", "#C70E7BFF", '#FC6882FF', "#C70E7BFF",
          '#FC6882FF', "#C70E7BFF", '#FC6882FF', 
          "#C70E7BFF", '#FC6882FF', "#C70E7BFF",
          '#A6E000FF', '#1BB6AFFF', "#C70E7BFF", '#6C6C9DFF', '#A6E000FF',
          '#1BB6AFFF', "#C70E7BFF", 
          '#A6E000FF','#1BB6AFFF', "#C70E7BFF", '#6C6C9DFF', 
          '#A6E000FF','#1BB6AFFF',  "#C70E7BFF", '#6C6C9DFF','#A6E000FF',
          '#1BB6AFFF', '#6C6C9DFF', "#C70E7BFF", 
          '#FC6882FF', "#C70E7BFF", '#6C6C9DFF',
          '#6C6C9DFF', "#C70E7BFF", 
          '#6C6C9DFF', 
          '#6C6C9DFF',  "#C70E7BFF",
          '#172869FF')

#colors time points
col1 <- c('#9E0142FF', '#9E0142FF', 
          '#D53E4FFF', '#D53E4FFF', '#D53E4FFF', '#D53E4FFF',
          '#F46D43FF', '#F46D43FF','#F46D43FF',
          '#F46D43FF','#F46D43FF','#F46D43FF',
          '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  
          '#FDAE61FF',  '#FDAE61FF',  
          '#FEE08BFF','#FEE08BFF','#FEE08BFF','#FEE08BFF',
          '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  
          '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF', 
          '#E6F598FF', '#E6F598FF', '#E6F598FF', 
          '#ABDDA4FF',  '#ABDDA4FF', 
          '#66C2A5FF',  
          '#3288BDFF',  '#3288BDFF',  
          '#5E4FA2FF')


DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, label=TRUE, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FNRMS/3_UMAP_projections.pdf", width=4, height=4, dpi=300)

DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, label=TRUE, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FNRMS/3_UMAP_projections_pop.pdf", width=4, height=4, dpi=300)

DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FNRMS/3_UMAP_projections_noleg.pdf", width=4, height=4, dpi=300)

DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
ggsave("./output/metadata/dev_muscle/SingleR_projection_FNRMS/3_UMAP_projections_pop_noleg.pdf", width=4, height=4, dpi=300)

DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012))
ggsave("./output/metadata/dev_muscle/SingleR_projection_FNRMS/3_UMAP_projections_legend.pdf", width=10, height=5, dpi=300)

DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col2, pt.size = 2, raster=TRUE, raster.dpi = c(1012, 1012))
ggsave("./output/metadata/dev_muscle/SingleR_projection_FNRMS/3_UMAP_projections_pop_legend.pdf", width=10, height=5, dpi=300)


#Bar plot clusters
# for each sample
plot_bar <-ggplot(FNRMS@meta.data, aes(x=name, fill= FNRMS@meta.data$SingleR.labels)) + 
  geom_bar(position = "fill", color="black") +
  scale_fill_manual(values = col1) +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_bar
ggsave("./output/metadata/dev_muscle/SingleR_projection_FNRMS/4_cell_composition.pdf", 
       width=10, height=6, dpi=300)

plot_bar <-ggplot(FNRMS@meta.data, aes(x=name, fill= FNRMS@meta.data$SingleR.labels)) + 
  geom_bar(position = "fill", color="black") +
  scale_fill_manual(values = col2) +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_bar
ggsave("./output/metadata/dev_muscle/SingleR_projection_FNRMS/4_cell_composition_pop.pdf", 
       width=10, height=6, dpi=300)


## Calculate number of cells within each cluster
library(data.table)
md <- FNRMS@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "./output/metadata/dev_muscle/SingleR_projection_FNRMS/5_Cluster_information_name.csv")



#save file
saveRDS(FNRMS, file = "./write/Danielli_Patel_Langenau_RPCA_ERMS_annotated_label.rds")
FNRMS <- readRDS("./write/Danielli_Patel_Langenau_RPCA_ERMS_annotated_label.rds")





