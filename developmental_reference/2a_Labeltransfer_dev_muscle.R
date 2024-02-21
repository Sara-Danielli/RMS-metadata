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
library(SCpubr)

#setwd('/Volumes/Sara_PhD/scRNAseq_data/')

######### load data ##########    
      FNRMS <- readRDS("../../Seurat_obj/Danielli_Patel_Langenau_RPCA_ERMS_20230713.rds")
      FPRMS.P3F1 <- readRDS("../../Seurat_obj/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_20230713.rds")
      FPRMS.P7F1 <- readRDS("../../Seurat_obj/Danielli_Patel_Langenau_RPCA_ARMS_P7F1_20230713.rds")
      muscle <- readRDS("processed reference muscle data.Rds")
      
      
      # reorganize labels to match Xi et al.
      Idents(muscle) = 'orig.ident.short'
      annot.ids <- c("Wk6.0_AP5" = "Wk5-6",
                     "Wk6.0_AP13" = "Wk5-6",
                     "Wk5.0_AP37" = "Wk5-6",
                     "Wk5.0_AP32" = "Wk5-6",
                     "Wk7.0_AP6" = "Wk6-7",
                     "Wk6.5_AP22" = "Wk6-7",
                     "Wk7.0_AP24" = "Wk6-7",
                     "Wk7.0_AP30" = "Wk6-7",
                     "Wk7.75_AP10" = "Wk7-8",
                     "Wk7.25_AP29" = "Wk7-8",
                     "Wk7.25_AP38" = "Wk7-8",
                     "Wk9_AP14" = "Wk9",
                     "Wk9_AP36" = "Wk9",
                     "Wk12_AP23" = "Wk12-14",
                     "Wk14_AP25" = "Wk12-14",
                     "Wk18_AP16" = "Wk17-18",
                     "Wk17_AP17" = "Wk17-18",
                     "Yr7_AP11" = "Postnatal",
                     "Yr11_AP34" = "Postnatal",
                     "Yr34_AP20" = "Postnatal",
                     "Yr42_AP28" = "Postnatal")
      muscle <- RenameIdents(muscle, annot.ids)
      muscle$annot.ids <- muscle@active.ident
     
      #define cell timepoint markers and select top 25 marker genes per time point
      muscle.markers <- FindAllMarkers(muscle, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3, test.use = 'negbinom')
      # muscle.markers.regressed <- FindAllMarkers(muscle, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3, test.use = 'negbinom',
      #                                            latent.vars = c("G2M.Score", "S.Score"))
      
      genes_per_cluster = 25
      top_genes <- muscle.markers %>%
        group_by(cluster) %>%
        slice_max(n = genes_per_cluster, order_by = avg_log2FC)
      
      top_gene_markers <- split(top_genes$gene, ceiling(seq_along(top_genes$gene) / genes_per_cluster))
      names(top_gene_markers) <- unique(top_genes$cluster)
      write.csv(x = data.frame(matrix(unlist(top_gene_markers),
                                      nrow = length(top_gene_markers), 
                                      byrow = T), 
                               row.names = names(top_gene_markers)), 
                file = "muscle timeline markers.csv")
      
      Idents(muscle) = 'orig.ident.short_cell.type'
      annot.ids <- c("Wk5.0_AP32_MP" = "Myogenic Progenitors",
                     "Wk5.0_AP37_MP" = "Myogenic Progenitors",
                     "Wk6.0_AP5_MP" = "Myogenic Progenitors",
                     "Wk6.0_AP13_MP" = "Myogenic Progenitors",
                     "Wk6.5_AP22_MB-MC" = "Myoblast-Myocytes",
                     "Wk6.5_AP22_MP" = "Myogenic Progenitors",
                     "Wk7.0_AP24_MB-MC" = "Myoblast-Myocytes",
                     "Wk7.0_AP24_MP" = "Myogenic Progenitors",
                     "Wk7.0_AP30_MB-MC" = "Myoblast-Myocytes",
                     "Wk7.0_AP30_MP" = "Myogenic Progenitors",
                     "Wk7.0_AP6_MB-MC" = "Myoblast-Myocytes",
                     "Wk7.0_AP6_MP" = "Myogenic Progenitors",
                     "Wk7.25_AP29_MB" = "Myoblasts",
                     "Wk7.25_AP29_MC" = "Myocytes",
                     "Wk7.25_AP29_MP" = "Myogenic Progenitors",
                     "Wk7.25_AP29_SkM.Mesen" = "Skeletal mesenchymal",
                     "Wk7.25_AP38_MB" = "Myoblasts",
                     "Wk7.25_AP38_MC" = "Myocytes",
                     "Wk7.25_AP38_MP" = "Myogenic Progenitors",
                     "Wk7.25_AP38_SkM.Mesen" = "Skeletal mesenchymal",
                     "Wk7.75_AP10_MB" = "Myoblasts",
                     "Wk7.75_AP10_MC" = "Myocytes",
                     "Wk7.75_AP10_MP" = "Myogenic Progenitors",
                     "Wk7.75_AP10_SkM.Mesen" = "Skeletal mesenchymal",
                     "Wk9_AP14_MB" = "Myoblasts",
                     "Wk9_AP14_MC" = "Myocytes",
                     "Wk9_AP14_MP" = "Myogenic Progenitors",
                     "Wk9_AP14_SkM.Mesen" = "Skeletal mesenchymal",
                     "Wk9_AP36_MB" = "Myoblasts",
                     "Wk9_AP36_MC" = "Myocytes",
                     "Wk9_AP36_MP" = "Myogenic Progenitors",
                     "Wk9_AP36_SkM.Mesen" = "Skeletal mesenchymal",
                     "Wk12_AP23_MB-MC" = "Myoblast-Myocytes",
                     "Wk12_AP23_MP" = "Myogenic Progenitors",
                     "Wk12_AP23_SkM.Mesen" = "Skeletal mesenchymal",
                     "Wk14_AP25_MB-MC" = "Myoblast-Myocytes",
                     "Wk14_AP25_SkM.Mesen" = "Skeletal mesenchymal",
                     "Wk14_AP25_MP" = "Myogenic Progenitors",
                     "Wk17_AP17_MB-MC" = "Myoblast-Myocytes",
                     "Wk17_AP17_MP" = "Myogenic Progenitors",
                     "Wk17_AP17_SkM.Mesen" = "Skeletal mesenchymal",
                     "Wk18_AP16_MB-MC" = "Myoblast-Myocytes",
                     "Wk18_AP16_MP" = "Myogenic Progenitors",
                     "Wk18_AP16_SkM.Mesen" = "Skeletal mesenchymal",
                     "Yr7_AP11_SC" = "Postnatal satellite cells",
                     "Yr11_AP34_SC" = "Postnatal satellite cells",
                     "Yr34_AP20_SC" = "Postnatal satellite cells",
                     "Yr42_AP28_SC" = "Postnatal satellite cells")
      muscle <- RenameIdents(muscle, annot.ids)
      muscle$cell.ids <- muscle@active.ident
      
      #define cell timepoint markers and select top 25 marker genes per time point
      muscle.cell.markers <- FindAllMarkers(muscle, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3, test.use = 'negbinom')
      # muscle.markers.regressed <- FindAllMarkers(muscle, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3, test.use = 'negbinom',
      #                                            latent.vars = c("G2M.Score", "S.Score"))
      
      genes_per_cluster = 25
      top_genes.cell <- muscle.cell.markers %>%
        group_by(cluster) %>%
        slice_max(n = genes_per_cluster, order_by = avg_log2FC)
      
      top_gene_markers.cell <- split(top_genes.cell$gene, ceiling(seq_along(top_genes.cell$gene) / genes_per_cluster))
      names(top_gene_markers.cell) <- unique(top_genes.cell$cluster)
      write.csv(x = data.frame(matrix(unlist(top_gene_markers.cell),
                                      nrow = length(top_gene_markers.cell), 
                                      byrow = T), 
                               row.names = names(top_gene_markers.cell)), 
                file = "muscle celltype markers.csv")
      
      # transform into single cell object for SingleR
      muscle.sce <- as.SingleCellExperiment(muscle)
    
      FNRMS@assays$ADT <- NULL
      FNRMS.sce <- as.SingleCellExperiment(FNRMS)
      
      FPRMS.P3F1@assays$ADT <- NULL
      FPRMS.P3F1.sce <- as.SingleCellExperiment(FPRMS.P3F1)
      
      FPRMS.P7F1@assays$ADT <- NULL
      FPRMS.P7F1.sce <- as.SingleCellExperiment(FPRMS.P7F1)

############################################################
      ## transfer muscle time labels onto RMS tumors
############################################################
      
# use SingleR to transfer labels
      pred.FNRMS <- SingleR(test=FNRMS.sce, ref=muscle.sce, labels=muscle.sce$annot.ids, genes = top_gene_markers)
      pred.FPRMS.P3F1 <- SingleR(test=FPRMS.P3F1.sce, ref=muscle.sce, labels=muscle.sce$annot.ids, genes = top_gene_markers)
      pred.FPRMS.P7F1 <- SingleR(test=FPRMS.P7F1.sce, ref=muscle.sce, labels=muscle.sce$annot.ids, genes = top_gene_markers)
      
      prop.table(table(pred.FNRMS$labels))*100
      prop.table(table(pred.FPRMS.P3F1$labels))*100
      prop.table(table(pred.FPRMS.P7F1$labels))*100
      
      write.csv(table(pred.FNRMS$labels), "0 fusion negative prediction.csv")
      write.csv(table(pred.FPRMS.P3F1$labels), "0 PAX3-FOXO1 prediction.csv")
      write.csv(table(pred.FPRMS.P7F1$labels), "0 PAX7-FOXO1 prediction.csv")
      
#plot heatmap
      pdf(file="./1_heatmap.pdf",
          width=9, height=5)
      print(plotScoreHeatmap(pred.FNRMS, order.by='label'))
      print(plotScoreHeatmap(pred.FPRMS.P3F1, order.by='label'))
      print(plotScoreHeatmap(pred.FPRMS.P7F1, order.by='label'))
      dev.off()
      
# plot  per-cell ???deltas??? = difference between the score for the assigned label and the median across all labels for each cell.
      # pdf(file="./2_deltas.pdf",
      #     width=6, height=10)
      # print(plotDeltaDistribution(pred.FPRMS, ncol = 3))
      # print(plotDeltaDistribution(pred.FPRMS2, ncol = 3))
      # print(plotDeltaDistribution(pred.FNRMS, ncol = 3))
      # dev.off()
      
# add predicted labels to Seurat object
      FNRMS[["SingleR.labels"]] <- pred.FNRMS$labels
      FPRMS.P3F1[["SingleR.labels"]] <- pred.FPRMS.P3F1$labels
      FPRMS.P7F1[["SingleR.labels"]] <- pred.FPRMS.P7F1$labels

# reorder labels
      FNRMS$SingleR.labels <- factor(x = FNRMS$SingleR.labels, 
                                     levels = c('Wk5-6',  'Wk6-7', 'Wk7-8', 'Wk9', 
                                                'Wk12-14', 'Wk17-18', 'Postnatal'))
      FPRMS.P3F1$SingleR.labels <- factor(x = FPRMS.P3F1$SingleR.labels, 
                                          levels = c('Wk5-6',  'Wk6-7', 'Wk7-8', 'Wk9', 
                                                     'Wk12-14', 'Wk17-18', 'Postnatal'))
      FPRMS.P7F1$SingleR.labels <- factor(x = FPRMS.P7F1$SingleR.labels , 
                                          levels = c('Wk5-6',  'Wk6-7', 'Wk7-8', 'Wk9', 
                                                     'Wk12-14', 'Wk17-18', 'Postnatal'))

#colors time points
col1 <- c('#9E0142FF', '#F46D43FF',  '#FDAE61FF', '#FFFFBFFF',  
          '#ABDDA4FF', '#66C2A5FF', '#3288BDFF', "#5E4FA2FF")

p1 <- DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + ggtitle('FN-RMS')
p2 <- DimPlot(FPRMS.P7F1, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + ggtitle('FP-RMS (PAX7::FOXO1)')
p3 <- DimPlot(FPRMS.P3F1, reduction = "umap_rpca", group.by = 'SingleR.labels', cols = col1, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + ggtitle('FP-RMS (PAX3::FOXO1')
p1| p2 | p3
ggsave("3_UMAP_projections.pdf", width=12.5, height=4, dpi=300)

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

plot_bar <-ggplot(FPRMS.P7F1@meta.data, aes(x=name, fill= FPRMS.P7F1@meta.data$SingleR.labels)) + 
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

plot_bar <-ggplot(FPRMS.P3F1@meta.data, aes(x=name, fill= FPRMS.P3F1@meta.data$SingleR.labels)) + 
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
 
## Calculate number of cells within each cluster
library(data.table)
md <- FNRMS@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_FN.csv")

md <- FPRMS.P3F1@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_P3F1.csv")

md <- FPRMS.P7F1@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_P7F1.csv")


############################################################
## transfer muscle cell labels onto RMS tumors
############################################################

# use SingleR to transfer labels
pred.cell.FNRMS <- SingleR(test=FNRMS.sce, ref=muscle.sce, labels=muscle.sce$cell.ids, genes = top_gene_markers.cell)
pred.cell.FPRMS.P3F1 <- SingleR(test=FPRMS.P3F1.sce, ref=muscle.sce, labels=muscle.sce$cell.ids, genes = top_gene_markers.cell)
pred.cell.FPRMS.P7F1 <- SingleR(test=FPRMS.P7F1.sce, ref=muscle.sce, labels=muscle.sce$cell.ids, genes = top_gene_markers.cell)

prop.table(table(pred.cell.FNRMS$labels))*100
prop.table(table(pred.cell.FPRMS.P3F1$labels))*100
prop.table(table(pred.cell.FPRMS.P7F1$labels))*100

write.csv(table(pred.cell.FNRMS$labels), "0 fusion negative cell prediction.csv")
write.csv(table(pred.cell.FPRMS.P3F1$labels), "0 PAX3-FOXO1 cell prediction.csv")
write.csv(table(pred.cell.FPRMS.P7F1$labels), "0 PAX7-FOXO1 cell prediction.csv")

write.csv(table(muscle$cell.ids, muscle$annot.ids), "0 time point cell props.csv")
#plot heatmap
pdf(file="./1_cell_heatmap.pdf",
    width=9, height=5)
print(plotScoreHeatmap(pred.cell.FNRMS, order.by='label'))
print(plotScoreHeatmap(pred.cell.FPRMS.P3F1, order.by='label'))
print(plotScoreHeatmap(pred.cell.FPRMS.P7F1, order.by='label'))
dev.off()

# plot  per-cell ???deltas??? = difference between the score for the assigned label and the median across all labels for each cell.
# pdf(file="./2_deltas.pdf",
#     width=6, height=10)
# print(plotDeltaDistribution(pred.cell.FPRMS, ncol = 3))
# print(plotDeltaDistribution(pred.cell.FPRMS2, ncol = 3))
# print(plotDeltaDistribution(pred.cell.FNRMS, ncol = 3))
# dev.off()

# add predicted labels to Seurat object
FNRMS[["SingleR.cell.labels"]] <- pred.cell.FNRMS$labels
FPRMS.P3F1[["SingleR.cell.labels"]] <- pred.cell.FPRMS.P3F1$labels
FPRMS.P7F1[["SingleR.cell.labels"]] <- pred.cell.FPRMS.P7F1$labels

# reorder labels
FNRMS$SingleR.cell.labels <- factor(x = FNRMS$SingleR.cell.labels, 
                               levels = c('Myogenic Progenitors',  'Myoblast-Myocytes', 'Myoblasts', 
                                          'Myocytes', 'Skeletal mesenchymal', 'Postnatal satellite cells'))
FPRMS.P3F1$SingleR.cell.labels <- factor(x = FPRMS.P3F1$SingleR.cell.labels, 
                                         levels = c('Myogenic Progenitors',  'Myoblast-Myocytes', 'Myoblasts', 
                                                    'Myocytes', 'Skeletal mesenchymal', 'Postnatal satellite cells'))
FPRMS.P7F1$SingleR.cell.labels <- factor(x = FPRMS.P7F1$SingleR.cell.labels , 
                                         levels = c('Myogenic Progenitors',  'Myoblast-Myocytes', 'Myoblasts', 
                                                    'Myocytes', 'Skeletal mesenchymal', 'Postnatal satellite cells'))
#colors time pointsdev
col2 <- c('#C7197CFF', '#F16882FF',  '#A4CF38FF', 
          '#1BB5AEFF', '#6D6C9DFF', '#262E67FF')

p1 <- DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col2, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + ggtitle('FN-RMS')
p2 <- DimPlot(FPRMS.P7F1, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col2, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + ggtitle('FP-RMS (PAX7::FOXO1)')
p3 <- DimPlot(FPRMS.P3F1, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col2, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + ggtitle('FP-RMS (PAX3::FOXO1')
p1| p2 | p3
ggsave("3_UMAP_cell_projections.pdf", width=12.5, height=4, dpi=300)

#Bar plot clusters
# for each sample
plot_bar <-ggplot(FNRMS@meta.data, aes(x=name, fill= FNRMS@meta.data$SingleR.cell.labels)) + 
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

plot_bar <-ggplot(FPRMS.P7F1@meta.data, aes(x=name, fill= FPRMS.P7F1@meta.data$SingleR.cell.labels)) + 
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

plot_bar <-ggplot(FPRMS.P3F1@meta.data, aes(x=name, fill= FPRMS.P3F1@meta.data$SingleR.cell.labels)) + 
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

## Calculate number of cells within each cluster
library(data.table)
md <- FNRMS@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.cell.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.cell.labels")] %>% dcast(., name ~ SingleR.cell.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_FN_cells.csv")

md <- FPRMS.P3F1@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.cell.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.cell.labels")] %>% dcast(., name ~ SingleR.cell.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_P3F1_cells.csv")

md <- FPRMS.P7F1@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.cell.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.cell.labels")] %>% dcast(., name ~ SingleR.cell.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_P7F1_cells.csv")


p1 <- DimPlot(muscle, reduction = "UMAP_no_regression", 
              group.by = 'annot.ids', cols = col1, label=F, pt.size = 3, 
              raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
p2 <- DimPlot(muscle, reduction = "UMAP_no_regression", 
              group.by = 'cell.ids', cols = col2, label=F, pt.size = 3, 
              raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend()
p1+p2
ggsave("S_UMAP_references.pdf", width=8.5, height=4, dpi=300)

#save file
# saveRDS(FPRMS.P3F1, file = "./Danielli_Patel_Langenau_RPCA_ARMS_P3F1_annotated_label.rds")
# saveRDS(FPRMS.P7F1, file = "./Danielli_Patel_Langenau_RPCA_ARMS_P7F1_annotated_label.rds")
# saveRDS(FNRMS, file = "./Danielli_Patel_Langenau_RPCA_ERMS_annotated_label.rds")


################################
FNRMS$SingleR.cell.labels <- factor(x = FNRMS$SingleR.cell.labels, 
                                    levels = c('Skeletal mesenchymal', 'Myogenic Progenitors',  'Myoblasts', 
                                               'Myoblast-Myocytes', 'Myocytes',  'Postnatal satellite cells'))
FPRMS.P3F1$SingleR.cell.labels <- factor(x = FPRMS.P3F1$SingleR.cell.labels, 
                                         levels = c('Skeletal mesenchymal', 'Myogenic Progenitors',  'Myoblasts', 
                                                    'Myoblast-Myocytes', 'Myocytes',  'Postnatal satellite cells'))
FPRMS.P7F1$SingleR.cell.labels <- factor(x = FPRMS.P7F1$SingleR.cell.labels , 
                                         levels = c('Skeletal mesenchymal', 'Myogenic Progenitors',  'Myoblasts', 
                                                    'Myoblast-Myocytes', 'Myocytes',  'Postnatal satellite cells'))

p1 <- SCpubr::do_AlluvialPlot(sample = FNRMS, 
                              first_group = "cluster_names_aggregate", 
                              last_group = "SingleR.cell.labels",
                              fill.by = "cluster_names_aggregate",
                              colors.use = c("Progenitor" = "#888F4B",
                                             "Proliferative" = "#FFAD72",
                                             "Differentiated" = "#8E0152",
                                             "IFN" = "#B496E6",
                                             "Ground" = "#FEE4CC"),
                              plot.title = "FN-RMS")

p2 <- SCpubr::do_AlluvialPlot(sample = FPRMS.P7F1, 
                              first_group = "cluster_names_aggregate", 
                              last_group = "SingleR.cell.labels",
                              fill.by = "cluster_names_aggregate",
                              colors.use = c("Progenitor" = "#888F4B",
                                             "Proliferative" = "#FFAD72",
                                             "Differentiated" = "#8E0152",
                                             "Neuronal" = "#B496E6",
                                             "Apoptosis" = "#D8D8E0",
                                             "Ground" = "#FEE4CC"),
                              plot.title = "FP-RMS (PAX7::FOXO1)")

p3 <- SCpubr::do_AlluvialPlot(sample = FPRMS.P3F1, 
                              first_group = "cluster_names_aggregate", 
                              last_group = "SingleR.cell.labels",
                              fill.by = "cluster_names_aggregate",
                              colors.use = c("Progenitor" = "#888F4B",
                                             "Proliferative" = "#FFAD72",
                                             "Differentiated" = "#8E0152",
                                             "Neuronal" = "#B496E6",
                                             "Apoptosis" = "#D8D8E0",
                                             "Ground" = "#FEE4CC"),
                              plot.title = "FP-RMS (PAX3::FOXO1)")

p1 + p2 + p3
ggplot2::ggsave("sankey plots.png", width = 21, height = 7, units = "in")

pdf("sankey plots.pdf", width = 21, height = 7)
print(p1+p2+p3)
dev.off()