rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratObject)
library(readxl)
library(data.table)
library(ggpubr)

# Set up environment ------------------------------------------------
base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
data_dir <- file.path(base_dir, 'data/DeMartino_Rds')
genelist_dir <- file.path(base_dir, 'list_final')

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/DeMartino')

plot_dir <- file.path(analysis_dir, 'plot')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
write_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(write_dir)){dir.create(write_dir, recursive = T)}

# color palette
col_subtype <- c('#95CECFFF', '#D3A2C2FF')
names(col_subtype) <- c('FP-RMS', 'FN-RMS')

# Load data ------------------------------------------------
# load RMS DeMartino
DeMartino <- readRDS(file.path(data_dir, "primary_filtered_malignant.rds"))

# rename subtypes
Idents(DeMartino) = "fusion_type"
new.cluster.ids.aggregate <- c('FP-RMS', 'FN-RMS', 'FP-RMS', 'FP-RMS')
names(new.cluster.ids.aggregate) <- levels(DeMartino)
DeMartino<- RenameIdents(DeMartino, new.cluster.ids.aggregate)
DeMartino[["subtype"]] <- Idents(object = DeMartino)

# Load RMS atlas signatures
Common_differentiated <- read_excel(file.path(genelist_dir, "Common_differentiated.xlsx"), col_names=FALSE)
Common_differentiated <- as.list(Common_differentiated)
Common_proliferative <- read_excel(file.path(genelist_dir,"Common_proliferative.xlsx"), col_names=FALSE)
Common_proliferative <- as.list(Common_proliferative)
Common_Stemcell <- read_excel(file.path(genelist_dir, "Common_Stemcell.xlsx"), col_names=FALSE)
Common_Stemcell <- as.list(Common_Stemcell)

signatures <- list(Common_differentiated, Common_proliferative, Common_Stemcell)
names(signatures) <- list('RMS_differentiated', 'RMS_proliferative', 'RMS_progenitor')

# Score data ------------------------------------------------
DeMartino <- AddModuleScore(object = DeMartino, features = Common_differentiated, name = 'Differentiated')
DeMartino <- AddModuleScore(object = DeMartino, features = Common_proliferative, name = 'Proliferative')
DeMartino <- AddModuleScore(object = DeMartino, features = Common_Stemcell, name = 'Progenitor')

DeMartino$DS.Difference <- DeMartino$Differentiated1 - DeMartino$Progenitor1


# Scaling data
DeMartino <- ScaleData(DeMartino)


# Plot data ------------------------------------------------

Vln_scores <- c("Differentiated1", "Proliferative1", "Progenitor1", "DS.Difference")
names(Vln_scores) <- c("Differentiated", "Proliferative", "Progenitor", "Muscle lineage score")

## Violin plots scores by sample
for (a in 1:length(Vln_scores)) {
  VlnPlot(DeMartino,
          features = Vln_scores[[a]], 
          group.by = 'tumor_id',  
          sort = 'decreasing',
          split.by = 'subtype',  pt.size=0,  
          cols = col_subtype) 
  ggsave(file.path(plot_dir, paste0("1_Vln_plot_scores_name_",names(Vln_scores)[a],".pdf")), width=6, height=3.5, dpi=300)
}

## Violin plots scores by subtype
for (a in 1:length(Vln_scores)) {
  VlnPlot(DeMartino,
          features = Vln_scores[[a]], 
          group.by = 'subtype',  
          sort = 'decreasing',
          pt.size=0,  
          cols = col_subtype) + 
    labs (y='Module score, AU', x='', title = names(Vln_scores)[a]) + 
    scale_fill_manual(values = col_subtype) + 
    geom_boxplot(outlier.shape=NA, width=0.1, fill="white") +
    stat_compare_means(aes(label = after_stat(p.signif)), 
                       method = 't.test', 
                       ref.group = 'FN-RMS',
                       size = 5, 
                       label.y.npc = 0.91, 
                       label.x.npc = 0.5) + 
    NoLegend() 
  ggsave(file.path(plot_dir, paste0("2_Vln_plot_scores_subtype_",names(Vln_scores)[a],".pdf")), width=3, height=3.5, dpi=300)
}
