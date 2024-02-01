rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratObject)
library(readxl)
library(data.table)
library(ggpubr)
library(cowplot)

# Set up environment ------------------------------------------------
base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
data_dir <- file.path(base_dir, 'data/DeMartino_Rds')
genelist_dir <- file.path(base_dir, 'list_final')

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/DeMartino')

resource_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources')
source(file.path(resource_dir, "Plot_style_v2.R"))


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
signatures <- read_excel(file.path(genelist_dir, "signatures_final.xlsx"), col_names=TRUE)
signatures <- as.list(signatures)
signatures <- sapply(signatures, function(x) x[!is.na(x)], simplify = FALSE)
names(signatures) <- list('Differentiated score', 'Proliferative score', 'Progenitor score')

# Score data ------------------------------------------------
DeMartino <- AddModuleScore(object = DeMartino, assay = 'RNA', features = signatures, name = names(signatures))

# rename metadata names of scores
# identify number of first column with metadata scores
col_start <- length(colnames(DeMartino@meta.data)) - length(names(signatures)) +1
# identify number of last column with metadata scores
col_end <- length(colnames(DeMartino@meta.data))
# rename columns with score name
colnames(DeMartino@meta.data)[col_start:col_end] <- names(signatures)

# Add Muscle lineage score
DeMartino$`Muscle lineage score` <- DeMartino$`Differentiated score` - DeMartino$`Progenitor score`

# Scaling data
DeMartino <- ScaleData(DeMartino)

# Save data ------------------------------------------------
saveRDS(DeMartino, file.path(base_dir, 'write/DeMartino_scored.rds'))


# Plot data ------------------------------------------------
scores <- c('Differentiated score', 'Proliferative score', 'Progenitor score', "Muscle lineage score")
titles <- c('Differentiated score', 'Proliferative score', 'Progenitor score', "Muscle lineage score")

p <- VlnPlot(DeMartino, features = scores, 
             group.by = 'tumor_id',  
             sort = 'decreasing',
             split.by = 'subtype', 
             cols = col_subtype,
             pt.size=0) 

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() +
    labs (y='Module score, AU', x='', title = titles[[i]]) + 
    scale_fill_manual(values = col_subtype) + 
    geom_boxplot(outlier.shape=NA, width=0.1, fill="white") + NoLegend() +
    theme_vln
}

plot <- p[[1]] + p[[2]] + p[[3]]+ p[[4]] 
plot

ggsave(file.path(plot_dir, paste0("1_Vln_plot_scores.pdf")), width=12, height=10)


## Violin plots scores by sample
for (a in 1:length(scores)) {
  VlnPlot(DeMartino,
          features = scores[[a]], 
          group.by = 'tumor_id',  
          sort = 'decreasing',
          split.by = 'subtype',  pt.size=0,  
          cols = col_subtype) + 
    labs (y='Module score, AU', x='', title = scores[a]) +  
    scale_fill_manual(values = col_subtype) + 
    geom_boxplot(outlier.shape=NA, width=0.1, fill="white") + 
    theme_vln
  ggsave(file.path(plot_dir, paste0("2_Vln_plot_scores_name_",scores[a],".pdf")), width=7, height=4, dpi=300)
}

## Violin plots scores by subtype
for (a in 1:length(scores)) {
  VlnPlot(DeMartino,
          features = scores[[a]], 
          group.by = 'subtype',  
          sort = 'decreasing',
          pt.size=0,  
          cols = col_subtype) + 
    labs (y='Module score, AU', x='', title = scores[a]) + 
    scale_fill_manual(values = col_subtype)  + 
    scale_fill_manual(values = col_subtype) + 
    geom_boxplot(outlier.shape=NA, width=0.1, fill="white") + NoLegend() +
    theme_vln +
    stat_compare_means(aes(label = after_stat(p.signif)), 
                       method = 't.test', 
                       #ref.group = 'FN-RMS',
                       size = 6, 
                       label.y.npc = 0.91, 
                       label.x.npc = 0.4) + 
    NoLegend() 
  ggsave(file.path(plot_dir, paste0("3_Vln_plot_scores_subtype_",scores[a],".pdf")), width=3, height=3.5, dpi=300)
}
