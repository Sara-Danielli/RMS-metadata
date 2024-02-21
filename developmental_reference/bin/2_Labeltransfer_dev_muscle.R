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
#library(scRNAseq)
library(SingleR)
library(SCpubr)

#options(Seurat.object.assay.version = 'v3')

# Set up environment ------------------------------------------------
base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')

output_dir <- file.path(base_dir, 'output/metadata/dev_muscle')
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

source(file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources/Plot_style_v2.R'))

# color palette
#colors populations
col_celltype <- c('#C70E7BFF', '#FC6882FF',  '#6C6C9DFF',  '#A6E000FF','#1BB6AFFF','#172869FF')
names(col_celltype) <- c('Myogenic progenitors', 'Myoblasts-myocyte', 'Skeletal mesenchymal', 'Myoblasts', 'Myocytes', 'Postnatal satellite cells')

#colors time points
col_timepoints <- c('#9E0142FF', '#D53E4FFF', '#F46D43FF', '#F46D43FF', '#FDAE61FF',  '#FEE08BFF',
                    '#FFFFBFFF',   '#E6F598FF', '#ABDDA4FF', '#66C2A5FF', '#3288BDFF', 
                    c(rep('#5E4FA2FF', 4)), 'grey')
names(col_timepoints) <- c("Wk5.0", "Wk6.0" ,  "Wk6.5", "Wk7.0", "Wk7.25",  "Wk7.75" , "Wk9" ,   "Wk12",   "Wk14", "Wk17" ,  "Wk18", 
                           "Yr7",    "Yr11",   "Yr34" ,  "Yr42", NA)


# Load seurat object data ------------------------------------------------
P3F1 <- readRDS(file.path(base_dir, "write/FPRMS_PAX3FOXO1_final_20240130.rds"))
P7F1 <- readRDS(file.path(base_dir, "write/FPRMS_PAX7FOXO1_final_20240130.rds"))
ERMS <- readRDS(file.path(base_dir, "write/FNRMS_final_20240130.rds"))

seurat_objects <- c(P3F1, P7F1, ERMS)
names(seurat_objects) <- c('P3F1', 'P7F1', 'FNRMS')

rm(P3F1, P7F1, ERMS)

muscle <- readRDS(file.path(base_dir, "write/Dev_muscle_Xietal.rds"))

# Load muscle signatures (Xi et al) ------------------------------------------------
muscle_signatures_celltype <- readRDS(file.path(output_dir, '7_markers_celltype_list.rds'))
  # select top 25 genes
  muscle_signatures_celltype <- lapply(muscle_signatures_celltype, head, 25)

muscle_signatures_timepoint <- readRDS(file.path(output_dir, '7_markers_timepoint_list.rds'))
  # select top 25 genes
  muscle_signatures_timepoint <- lapply(muscle_signatures_timepoint, head, 25)

# Transform objects into single cell  ------------------------------------------------      
muscle.sce <- as.SingleCellExperiment(muscle)

sc_objects <- list()

for (i in seq_along(seurat_objects)) {      
sc_objects[[i]] <- as.SingleCellExperiment(seurat_objects[[i]])
}
names(sc_objects) <- names(seurat_objects)


# Transfer muscle time labels onto RMS tumors ------------------------------------------------      
pred.labels <- list()
prop.table <- list()

for (i in seq_along(sc_objects)) {
# use SingleR to transfer labels
pred.labels[[i]] <- SingleR(test=sc_objects[[i]], ref=muscle.sce, labels=muscle.sce$celltype, genes = muscle_signatures_celltype)
#prop.table[[i]] <- table(pred.labels[[i]]$labels) %>% mutate(freq = n/)
#write.csv(prop.table[[i]], file.path(output_dir, paste0("0_", names(sc_objects)[i], "_prediction.csv")))

pdf(file.path(output_dir, paste0("1_", names(sc_objects)[i], "_score_heatmap_celltype.pdf")), width=9, height=5)
print(plotScoreHeatmap(pred.labels[[i]], order.by='label'))
dev.off()
}


for (i in seq_along(sc_objects)) {
  # use SingleR to transfer labels
  pred.labels[[i]] <- SingleR(test=sc_objects[[i]], ref=muscle.sce, labels=muscle.sce$timepoint, genes = muscle_signatures_timepoint)
  #prop.table[[i]] <- table(pred.labels[[i]]$labels) %>% mutate(freq = n/)
  #write.csv(prop.table[[i]], file.path(output_dir, paste0("0_", names(sc_objects)[i], "_prediction.csv")))
  
  pdf(file.path(output_dir, paste0("1_", names(sc_objects)[i], "_score_heatmap.pdf")), width=9, height=5)
  print(plotScoreHeatmap(pred.labels[[i]], order.by='label'))
  dev.off()
}
      

# add predicted labels to Seurat object
seurat_objects[[i]][["SingleR.labels"]] <- pred.labels[[i]]$labels

# reorder labels
seurat_objects[[i]]$SingleR.labels <- factor(x = seurat_objects[[i]]$SingleR.labels, 
                                             levels = c("Wk5.0", "Wk6.0" ,  "Wk6.5", "Wk7.0", "Wk7.25",  "Wk7.75" , "Wk9" ,   "Wk12",   "Wk14", "Wk17" ,  "Wk18", 
                                                        "Yr7",    "Yr11",   "Yr34" ,  "Yr42"))
# plot
DimPlot(seurat_objects[[i]], reduction = "umap_rpca", group.by = 'SingleR.labels', 
        cols = col1, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + 
  NoLegend() + ggtitle(names(seurat_objects)[i])
ggsave(file.path(output_dir, paste0("3_UMAP_projections", names(seurat_objects)[i], ".pdf", width=4, height=4)))

#Bar plot clusters
# for each sample
plot_bar <-ggplot(ERMS@meta.data, aes(x=name, fill= ERMS@meta.data$SingleR.labels)) + 
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

plot_bar <-ggplot(P7F1@meta.data, aes(x=name, fill= P7F1@meta.data$SingleR.labels)) + 
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

plot_bar <-ggplot(P3F1@meta.data, aes(x=name, fill= P3F1@meta.data$SingleR.labels)) + 
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
md <- ERMS@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_FN.csv")

md <- P3F1@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_P3F1.csv")

md <- P7F1@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.labels")] %>% dcast(., name ~ SingleR.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_P7F1.csv")

}



# Transfer muscle cell type labels onto RMS tumors ------------------------------------------------    

# use SingleR to transfer labels
pred.cell.ERMS <- SingleR(test=ERMS.sce, ref=muscle.sce, labels=muscle.sce$cell.ids, genes = muscle_signatures_celltype)
pred.cell.P3F1 <- SingleR(test=P3F1.sce, ref=muscle.sce, labels=muscle.sce$cell.ids, genes = muscle_signatures_celltype)
pred.cell.P7F1 <- SingleR(test=P7F1.sce, ref=muscle.sce, labels=muscle.sce$cell.ids, genes = muscle_signatures_celltype)

prop.table(table(pred.cell.ERMS$labels))*100
prop.table(table(pred.cell.P3F1$labels))*100
prop.table(table(pred.cell.P7F1$labels))*100

write.csv(table(pred.cell.ERMS$labels), "0 fusion negative cell prediction.csv")
write.csv(table(pred.cell.P3F1$labels), "0 PAX3-FOXO1 cell prediction.csv")
write.csv(table(pred.cell.P7F1$labels), "0 PAX7-FOXO1 cell prediction.csv")

write.csv(table(muscle$cell.ids, muscle$annot.ids), "0 time point cell props.csv")
#plot heatmap
pdf(file="./1_cell_heatmap.pdf",
    width=9, height=5)
print(plotScoreHeatmap(pred.cell.ERMS, order.by='label'))
print(plotScoreHeatmap(pred.cell.P3F1, order.by='label'))
print(plotScoreHeatmap(pred.cell.P7F1, order.by='label'))
dev.off()

# plot  per-cell ???deltas??? = difference between the score for the assigned label and the median across all labels for each cell.
# pdf(file="./2_deltas.pdf",
#     width=6, height=10)
# print(plotDeltaDistribution(pred.cell.FPRMS, ncol = 3))
# print(plotDeltaDistribution(pred.cell.FPRMS2, ncol = 3))
# print(plotDeltaDistribution(pred.cell.ERMS, ncol = 3))
# dev.off()

# add predicted labels to Seurat object
ERMS[["SingleR.cell.labels"]] <- pred.cell.ERMS$labels
P3F1[["SingleR.cell.labels"]] <- pred.cell.P3F1$labels
P7F1[["SingleR.cell.labels"]] <- pred.cell.P7F1$labels

# reorder labels
ERMS$SingleR.cell.labels <- factor(x = ERMS$SingleR.cell.labels, 
                               levels = c('Myogenic Progenitors',  'Myoblast-Myocytes', 'Myoblasts', 
                                          'Myocytes', 'Skeletal mesenchymal', 'Postnatal satellite cells'))
P3F1$SingleR.cell.labels <- factor(x = P3F1$SingleR.cell.labels, 
                                         levels = c('Myogenic Progenitors',  'Myoblast-Myocytes', 'Myoblasts', 
                                                    'Myocytes', 'Skeletal mesenchymal', 'Postnatal satellite cells'))
P7F1$SingleR.cell.labels <- factor(x = P7F1$SingleR.cell.labels , 
                                         levels = c('Myogenic Progenitors',  'Myoblast-Myocytes', 'Myoblasts', 
                                                    'Myocytes', 'Skeletal mesenchymal', 'Postnatal satellite cells'))
#colors time pointsdev
col2 <- c('#C7197CFF', '#F16882FF',  '#A4CF38FF', 
          '#1BB5AEFF', '#6D6C9DFF', '#262E67FF')

p1 <- DimPlot(ERMS, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col2, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + ggtitle('FN-RMS')
p2 <- DimPlot(P7F1, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col2, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + ggtitle('FP-RMS (PAX7::FOXO1)')
p3 <- DimPlot(P3F1, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col2, label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + ggtitle('FP-RMS (PAX3::FOXO1')
p1| p2 | p3
ggsave("3_UMAP_cell_projections.pdf", width=12.5, height=4, dpi=300)

#Bar plot clusters
# for each sample
plot_bar <-ggplot(ERMS@meta.data, aes(x=name, fill= ERMS@meta.data$SingleR.cell.labels)) + 
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

plot_bar <-ggplot(P7F1@meta.data, aes(x=name, fill= P7F1@meta.data$SingleR.cell.labels)) + 
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

plot_bar <-ggplot(P3F1@meta.data, aes(x=name, fill= P3F1@meta.data$SingleR.cell.labels)) + 
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
md <- ERMS@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.cell.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.cell.labels")] %>% dcast(., name ~ SingleR.cell.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_FN_cells.csv")

md <- P3F1@meta.data %>% as.data.table
md[, .N, by = c("name", "SingleR.cell.labels")]
table_PDX <- md[, .N, by = c("name", "SingleR.cell.labels")] %>% dcast(., name ~ SingleR.cell.labels, value.var = "N")
write.csv(table_PDX, "5_Cluster_information_name_P3F1_cells.csv")

md <- P7F1@meta.data %>% as.data.table
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
# saveRDS(P3F1, file = "./Danielli_Patel_Langenau_RPCA_ARMS_P3F1_annotated_label.rds")
# saveRDS(P7F1, file = "./Danielli_Patel_Langenau_RPCA_ARMS_P7F1_annotated_label.rds")
# saveRDS(ERMS, file = "./Danielli_Patel_Langenau_RPCA_ERMS_annotated_label.rds")


################################
ERMS$SingleR.cell.labels <- factor(x = ERMS$SingleR.cell.labels, 
                                    levels = c('Skeletal mesenchymal', 'Myogenic Progenitors',  'Myoblasts', 
                                               'Myoblast-Myocytes', 'Myocytes',  'Postnatal satellite cells'))
P3F1$SingleR.cell.labels <- factor(x = P3F1$SingleR.cell.labels, 
                                         levels = c('Skeletal mesenchymal', 'Myogenic Progenitors',  'Myoblasts', 
                                                    'Myoblast-Myocytes', 'Myocytes',  'Postnatal satellite cells'))
P7F1$SingleR.cell.labels <- factor(x = P7F1$SingleR.cell.labels , 
                                         levels = c('Skeletal mesenchymal', 'Myogenic Progenitors',  'Myoblasts', 
                                                    'Myoblast-Myocytes', 'Myocytes',  'Postnatal satellite cells'))

p1 <- SCpubr::do_AlluvialPlot(sample = ERMS, 
                              first_group = "cluster_names_aggregate", 
                              last_group = "SingleR.cell.labels",
                              fill.by = "cluster_names_aggregate",
                              colors.use = c("Progenitor" = "#888F4B",
                                             "Proliferative" = "#FFAD72",
                                             "Differentiated" = "#8E0152",
                                             "IFN" = "#B496E6",
                                             "Ground" = "#FEE4CC"),
                              plot.title = "FN-RMS")

p2 <- SCpubr::do_AlluvialPlot(sample = P7F1, 
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

p3 <- SCpubr::do_AlluvialPlot(sample = P3F1, 
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