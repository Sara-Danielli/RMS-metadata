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
library(stringr)

setwd('/Volumes/Sara_PhD/scRNAseq_data/RMS/Xi_Cell_stem_cell_2020')

##################################################################
#  (1) Load human muscle dataset Xi et al. 2020:  #
##################################################################

# load human development datasets Xi et al. 2020 Stem Cell
week5_6 <- read.table(file = './myogenic_datasets/week5_6/exprMatrix.tsv', sep = '\t', header=TRUE, row.names=1)
week5_6_meta <- read.table(file = './myogenic_datasets/week5_6/meta.tsv', sep = '\t', header=TRUE, row.names=1)
week5_6 <- CreateSeuratObject(counts = week5_6, project = "Xi_2020")
week5_6 <- AddMetaData(week5_6, metadata = week5_6_meta)

week6_7 <- read.table(file = './myogenic_datasets/week6_7/exprMatrix.tsv', sep = '\t', header=TRUE, row.names=1)
week6_7_meta <- read.table(file = './myogenic_datasets/week6_7/meta.tsv', sep = '\t', header=TRUE, row.names=1)
week6_7 <- CreateSeuratObject(counts = week6_7, project = "Xi_2020")
week6_7 <- AddMetaData(week6_7, metadata = week6_7_meta)

week7_8 <- read.table(file = './myogenic_datasets/week7_8/exprMatrix.tsv', sep = '\t', header=TRUE, row.names=1)
week7_8_meta <- read.table(file = './myogenic_datasets/week7_8/meta.tsv', sep = '\t', header=TRUE, row.names=1)
week7_8 <- CreateSeuratObject(counts = week7_8, project = "Xi_2020")
week7_8 <- AddMetaData(week7_8, metadata = week7_8_meta)

week9 <- read.table(file = './myogenic_datasets/week9/exprMatrix.tsv', sep = '\t', header=TRUE, row.names=1)
week9_meta<- read.table(file = './myogenic_datasets/week9/meta.tsv', sep = '\t', header=TRUE, row.names=1)
week9 <- CreateSeuratObject(counts = week9, project = "Xi_2020")
week9 <- AddMetaData(week9, metadata = week9_meta)

week12_14 <- read.table(file = './myogenic_datasets/week12_14/exprMatrix.tsv', sep = '\t', header=TRUE, row.names=1)
week12_14_meta<- read.table(file = './myogenic_datasets/week12_14/meta.tsv', sep = '\t', header=TRUE, row.names=1)
week12_14 <- CreateSeuratObject(counts = week12_14, project = "Xi_2020")
week12_14 <- AddMetaData(week12_14, metadata = week12_14_meta)

week17_18 <- read.table(file = './myogenic_datasets/week17_18/exprMatrix.tsv', sep = '\t', header=TRUE, row.names=1)
week17_18_meta<- read.table(file = './myogenic_datasets/week17_18/meta.tsv', sep = '\t', header=TRUE, row.names=1)
week17_18 <- CreateSeuratObject(counts = week17_18, project = "Xi_2020")
week17_18 <- AddMetaData(week17_18, metadata = week17_18_meta)

juvenile <- read.table(file = './myogenic_datasets/juvenile/exprMatrix.tsv', sep = '\t', header=TRUE, row.names=1)
juvenile_meta <- read.table(file = './myogenic_datasets/juvenile/meta.tsv', sep = '\t', header=TRUE, row.names=1)
juvenile <- CreateSeuratObject(counts = juvenile, project = "Xi_2020")
juvenile  <- AddMetaData(juvenile , metadata = juvenile_meta)

adult <- read.table(file = './myogenic_datasets/adult/exprMatrix.tsv', sep = '\t', header=TRUE, row.names=1)
adult_meta <- read.table(file = './myogenic_datasets/adult/meta.tsv', sep = '\t', header=TRUE, row.names=1)
adult <- CreateSeuratObject(counts = adult, project = "Xi_2020")
adult <- AddMetaData(adult, metadata = adult_meta)



merged_Xietal <- merge(week5_6, y = list(week6_7, week7_8, week9, week12_14, week17_18, juvenile, adult),
                       add.cell.ids = c('week5_6', 'week6_7', 'week7_8', 'week9', 'week12_14', 'week17_18', 'juvenile', 'adult'),
                       project = "development_Xi_2020")

Idents(merged_Xietal) <- 'orig.ident.short'

DefaultAssay(merged_Xietal) <- "RNA"

# Define mitochondrial genes
merged_Xietal[["percent.mt"]] <- PercentageFeatureSet(merged_Xietal, pattern = "^MT-")

#metadata
  # add cell type info
  cellname <- colnames(merged_Xietal)
  celltype <- unlist(str_split(merged_Xietal@meta.data$Cell.Type, "-", n = 1))
  #create dataframe
  df <- data.frame(cellname, celltype)
  # add info as metadata  
  AddMetaData(merged_Xietal, metadata = df)
  
#save object  
saveRDS(merged_Xietal, file = "/Volumes/Sara_PhD/scRNAseq_data/write/Dev_muscle_Xietal.rds")
  
merged_Xietal <- readRDS("/Volumes/Sara_PhD/scRNAseq_data/write/Dev_muscle_Xietal.rds")

##################################################################
#  (2) Downstream analyses:  #
##################################################################

## Data normalization
merged_Xietal <- NormalizeData(merged_Xietal, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable features
merged_Xietal <- FindVariableFeatures(merged_Xietal, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_Xietal), 10)
# Scaling data
merged_Xietal <- ScaleData(merged_Xietal, features = VariableFeatures(merged_Xietal))

## Dimension reduction
merged_Xietal
merged_Xietal <- RunPCA(merged_Xietal, reduction.name = "pca_no_regression", reduction.key = "pca_no_regression")

# reorder labels
merged_Xietal$orig.ident.short_cell.type <- factor(x = merged_Xietal$orig.ident.short_cell.type, 
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

#colors time points
col1 <- c('#9E0142FF', '#9E0142FF', 
          '#D53E4FFF', '#D53E4FFF', '#D53E4FFF', '#D53E4FFF',
          '#F46D43FF', '#F46D43FF','#F46D43FF',
          '#F46D43FF','#F46D43FF','#F46D43FF',
          '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF',  
          '#FDAE61FF',  '#FDAE61FF',  '#FDAE61FF', 
          '#FEE08BFF','#FEE08BFF','#FEE08BFF','#FEE08BFF',
          '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF',  
          '#FFFFBFFF',  '#FFFFBFFF',  '#FFFFBFFF', 
          '#E6F598FF', '#E6F598FF', '#E6F598FF', 
          '#ABDDA4FF',  '#ABDDA4FF', '#ABDDA4FF', 
          '#66C2A5FF',  '#66C2A5FF',  '#66C2A5FF', 
          '#3288BDFF',   '#3288BDFF',  '#3288BDFF',
          '#5E4FA2FF',  '#5E4FA2FF', '#5E4FA2FF', '#5E4FA2FF')


#colors populations
col2 <- c('#C70E7BFF', '#C70E7BFF', 
          '#C70E7BFF', '#C70E7BFF', '#FC6882FF', '#C70E7BFF',
          '#FC6882FF', '#C70E7BFF', '#FC6882FF', 
          '#C70E7BFF', '#FC6882FF', '#C70E7BFF',
          '#A6E000FF', '#1BB6AFFF', '#C70E7BFF', '#6C6C9DFF', '#A6E000FF',
          '#1BB6AFFF', '#C70E7BFF', '#6C6C9DFF',
          '#A6E000FF','#1BB6AFFF', '#C70E7BFF', '#6C6C9DFF', 
          '#A6E000FF','#1BB6AFFF',  '#C70E7BFF', '#6C6C9DFF','#A6E000FF',
          '#1BB6AFFF', '#6C6C9DFF', '#C70E7BFF', 
          '#FC6882FF', '#C70E7BFF', '#6C6C9DFF',
          '#FC6882FF', '#6C6C9DFF', '#C70E7BFF', 
          '#FC6882FF', '#C70E7BFF', '#6C6C9DFF', 
          '#FC6882FF', '#C70E7BFF', '#6C6C9DFF', 
          '#172869FF', '#172869FF', '#172869FF', '#172869FF')


DimPlot(merged_Xietal, reduction = "pca_no_regression", group.by = 'orig.ident.short_cell.type', cols = col1) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/1_PCA_celltype.pdf", 
       width=10, height=5, dpi=300)


ElbowPlot(merged_Xietal, reduction='pca_no_regression')
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/2_Elbowplot.pdf", 
       width=6, height=4, dpi=300)


## Cluster cells
merged_Xietal <- FindNeighbors(merged_Xietal, reduction = 'pca_no_regression', dims = 1:10)

# UMAP
merged_Xietal <- RunUMAP(merged_Xietal, dims = 1:10, reduction = 'pca_no_regression', reduction.name = 'UMAP_no_regression')

## Plots
DimPlot(merged_Xietal, reduction = "UMAP_no_regression", group.by = 'orig.ident.short_cell.type', cols = col1) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/3_UMAP_celltype_col1.pdf", 
       width=12, height=5, dpi=300)

DimPlot(merged_Xietal, reduction = "UMAP_no_regression", group.by = 'orig.ident.short_cell.type', cols=col2) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/3_UMAP_celltype_col2.pdf", 
       width=12, height=5, dpi=300)


plot_bar <-ggplot(merged_Xietal@meta.data, aes(x=orig.ident.short, fill= (merged_Xietal@meta.data$orig.ident.short_cell.type))) + 
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
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/4_cell_composition_pop.pdf", 
       width=10, height=4, dpi=300)

## Calculate number of cells within each cluster
library(data.table)
md <- merged_Xietal@meta.data %>% as.data.table
md[, .N, by = c("orig.ident.short", "orig.ident.short_cell.type")]
table_PDX <- md[, .N, by = c("orig.ident.short", "orig.ident.short_cell.type")] %>% dcast(., orig.ident.short ~ orig.ident.short_cell.type, value.var = "N")
write.csv(table_PDX, "/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/5_Cluster_information_name.csv")






