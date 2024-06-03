rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(cowplot)
library(ggpubr)
library(tidyverse)

# Set up environment ------------------------------------------------
base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
data_dir <- file.path(base_dir, 'RMS/Patel_2022/FFPE_RNAseq')
plot_dir <- file.path(base_dir, 'output/metadata/Pseudobulk')

if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

source(file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources/Plot_style_v2.R'))


# Load FFPE RNAseq datasets ------------------------------------------------

## ! Excel matrix contained some numbers that were not integers --> approximated

# load TPM data
dataset <- read.delim(file.path(data_dir, "RMS13_FFPE_RNA-seq_tpm_data.txt"))

rownames(dataset) <- make.names(dataset$geneSymbol, unique = TRUE)

# set gene names as row header
dataset2 <- dataset[, -1]
dataset3 <- dataset2[, -1]
dataset4 <- dataset3[, -1]
dataset5 <- dataset4[, -1]

# load metadata 
FFPE_meta <- read_excel(file.path(data_dir,"metadata.xlsx"), col_names=TRUE)
  # set first column as row header
FFPE_meta2 <- FFPE_meta %>% remove_rownames %>% column_to_rownames(var="patient_id")

# Create Seurat object
FFPE_seurat <- CreateSeuratObject(counts = dataset5)
FFPE_seurat <- AddMetaData(FFPE_seurat, metadata = FFPE_meta2)

# log-transform (as data are already TPM, do not perform "NormalizeData")
counts <- GetAssayData(FFPE_seurat, slot = "counts")
counts_log <-log1p(counts)

# Re-create Seurat object with log-transformed counts
FFPE_seurat2 <- CreateSeuratObject(counts = counts_log)
FFPE_seurat[["RNA"]]$data <- FFPE_seurat2[["RNA"]]$counts

# Reoder levels
FFPE_seurat$treatment <- factor(x = FFPE_seurat$treatment, levels = c('diagnostic biopsy', 'delayed resection'))


# Score datasets with gene signatures  ------------------------------------------------

# Load common signatures
gene_list <- as.list(read_excel(file.path(base_dir, "list_final/RMS_atlas_gene_lists.xlsx"), col_names=TRUE))
  # remove NAs
  df_clean <- lapply(gene_list, na.omit)
  df_list_clean <- as.list(df_clean)
  remove_na_action <- function(x) {
    attributes(x)$na.action <- NULL
    x
  }
  gene_list <- lapply(df_list_clean, function(x) remove_na_action(na.omit(x)))


# Add scores using "Addmodulescore" fct
  FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = gene_list, name = names(gene_list))
  
  # rename metadata names of scores
  col_start <- length(colnames(FFPE_seurat@meta.data)) - length(names(gene_list)) + 1
  # identify number of last column with metadata scores
  col_end <- length(colnames(FFPE_seurat@meta.data))
  # rename columns with score name
  colnames(FFPE_seurat@meta.data)[col_start:col_end] <- names(gene_list)
  
  FFPE_seurat$Muscle_lineage_score <- FFPE_seurat$Differentiated- FFPE_seurat$Progenitor
  
  
  # Scale data
ScaleData(FFPE_seurat)

# save object
saveRDS(FFPE_seurat, file.path(base_dir, 'write/FFPE_patients_processed.rds'))

      
# Violin plots ------------------------------------------------
      compare_means(Proliferative ~ treatment,  data = FFPE_seurat[[]], paired = TRUE)
      p2 <- VlnPlot(FFPE_seurat,
              features = 'Proliferative',
              group.by = 'treatment',
              pt.size=0,
              cols = c( 'gray96', 'gray96', 'lightseagreen', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat@meta.data$subtype,group=FFPE_seurat@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.signif)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      ggsave(file.path(plot_dir, "1_Vln_plot_prolifer_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      
      
      compare_means(Differentiated ~ treatment,  data = FFPE_seurat[[]], paired = TRUE)
      p3 <- VlnPlot(FFPE_seurat,
              features = 'Differentiated',
              group.by = 'treatment',
              pt.size=0,
              cols = c( 'gray96', 'gray96', 'lightseagreen', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat@meta.data$subtype, group=FFPE_seurat@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.signif)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      ggsave(file.path(plot_dir, "2_Vln_plot_differentiated_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      
      
      
      compare_means(Progenitor ~ treatment,  data = FFPE_seurat[[]])
      p1 <- VlnPlot(FFPE_seurat,
              features = 'Progenitor',
              group.by = 'treatment',
              pt.size=0,
              cols = c( 'gray96', 'gray96', 'lightseagreen', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat@meta.data$subtype, group=FFPE_seurat@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.signif)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      ggsave(file.path(plot_dir, "3_Vln_plot_progenitor_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      
      
      compare_means(Muscle_lineage_score ~ treatment,  data = FFPE_seurat[[]])
      VlnPlot(FFPE_seurat,
              features = 'Muscle_lineage_score',
              group.by = 'treatment',
              pt.size=0,
              cols = c( 'gray96', 'gray96', 'lightseagreen', 'maroon'),
      ) + theme_vln +
        geom_point(aes(fill=FFPE_seurat@meta.data$subtype, group=FFPE_seurat@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.signif)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()
      ggsave(file.path(plot_dir, "4_Vln_plot_muscle-score_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      plot_grid(p1, p2, p3, ncol = 3, align = "h", axis = "tb")
      ggsave(file.path(plot_dir,"5_VlnPlot_combined.pdf"), width=7, height=4)
      
 
      
# Violin plots FN-RMS only ------------------------------------------------      
      # subset
      FFPE_seurat_ERMS <- subset(FFPE_seurat, subset = subtype == 'ERMS')
      # Scale data
      ScaleData(FFPE_seurat_ERMS)
      
      
      # Violin plots ------------------------------------------------
      compare_means(Proliferative ~ treatment,  data = FFPE_seurat_ERMS[[]], paired = TRUE)
      p2 <- VlnPlot(FFPE_seurat_ERMS,
                    features = 'Proliferative',
                    group.by = 'treatment',
                    pt.size=0,
                    cols = c( 'gray96', 'gray96', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat_ERMS@meta.data$subtype,group=FFPE_seurat_ERMS@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat_ERMS@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.format)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      
      
      compare_means(Differentiated ~ treatment,  data = FFPE_seurat_ERMS[[]], paired = TRUE)
      p3 <- VlnPlot(FFPE_seurat_ERMS,
                    features = 'Differentiated',
                    group.by = 'treatment',
                    pt.size=0,
                    cols = c( 'gray96', 'gray96', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat_ERMS@meta.data$subtype, group=FFPE_seurat_ERMS@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat_ERMS@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.format)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      
      
      
      compare_means(Progenitor ~ treatment,  data = FFPE_seurat_ERMS[[]])
      p1 <- VlnPlot(FFPE_seurat_ERMS,
                    features = 'Progenitor',
                    group.by = 'treatment',
                    pt.size=0,
                    cols = c( 'gray96', 'gray96', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat_ERMS@meta.data$subtype, group=FFPE_seurat_ERMS@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat_ERMS@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.format)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      
      
      compare_means(Muscle_lineage_score ~ treatment,  data = FFPE_seurat_ERMS[[]])
      VlnPlot(FFPE_seurat_ERMS,
              features = 'Muscle_lineage_score',
              group.by = 'treatment',
              pt.size=0,
              cols = c( 'gray96', 'gray96', 'maroon'),
      ) + theme_vln +
        geom_point(aes(fill=FFPE_seurat_ERMS@meta.data$subtype, group=FFPE_seurat_ERMS@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat_ERMS@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.format)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()
      ggsave(file.path(plot_dir, "6_FNRMS__Vln_plot_muscle-score_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      plot_grid(p1, p2, p3, ncol = 3, align = "h", axis = "tb")
      ggsave(file.path(plot_dir,"7_FNRMS_VlnPlot.pdf"), width=7, height=4.5)
      
      write.csv(FFPE_seurat_ERMS@meta.data, file.path(plot_dir, '7_FNRMS_raw_data.csv'))
      
      
      
# Violin plots FP-RMS only ------------------------------------------------      
      # subset
      FFPE_seurat_ARMS <- subset(FFPE_seurat, subset = subtype == 'ARMS')
      # Scale data
      ScaleData(FFPE_seurat_ARMS)
      
      
      # Violin plots ------------------------------------------------
      compare_means(Proliferative ~ treatment,  data = FFPE_seurat_ARMS[[]], paired = TRUE)
      p2 <- VlnPlot(FFPE_seurat_ARMS,
                    features = 'Proliferative',
                    group.by = 'treatment',
                    pt.size=0,
                    cols = c( 'gray96', 'gray96', 'lightseagreen'),
      ) + 
        geom_point(aes(fill=FFPE_seurat_ARMS@meta.data$subtype,group=FFPE_seurat_ARMS@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat_ARMS@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.signif)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      
      
      compare_means(Differentiated ~ treatment,  data = FFPE_seurat_ARMS[[]], paired = TRUE)
      p3 <- VlnPlot(FFPE_seurat_ARMS,
                    features = 'Differentiated',
                    group.by = 'treatment',
                    pt.size=0,
                    cols = c( 'gray96', 'gray96', 'lightseagreen'),
      ) + 
        geom_point(aes(fill=FFPE_seurat_ARMS@meta.data$subtype, group=FFPE_seurat_ARMS@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat_ARMS@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.signif)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      
      
      
      compare_means(Progenitor ~ treatment,  data = FFPE_seurat_ARMS[[]])
      p1 <- VlnPlot(FFPE_seurat_ARMS,
                    features = 'Progenitor',
                    group.by = 'treatment',
                    pt.size=0,
                    cols = c( 'gray96', 'gray96', 'lightseagreen'),
      ) + 
        geom_point(aes(fill=FFPE_seurat_ARMS@meta.data$subtype, group=FFPE_seurat_ARMS@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat_ARMS@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.signif)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()      
      
      
      compare_means(Muscle_lineage_score ~ treatment,  data = FFPE_seurat_ARMS[[]])
      VlnPlot(FFPE_seurat_ARMS,
              features = 'Muscle_lineage_score',
              group.by = 'treatment',
              pt.size=0,
              cols = c( 'gray96', 'gray96', 'lightseagreen'),
      ) + theme_vln +
        geom_point(aes(fill=FFPE_seurat_ARMS@meta.data$subtype, group=FFPE_seurat_ARMS@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat_ARMS@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(aes(label = after_stat(p.signif)), 
                           method = 't.test', 
                           size = 5, 
                           label.y.npc = 0.95, 
                           label.x.npc = 0.5, 
        ) + NoLegend()
      ggsave(file.path(plot_dir, "6_FNRMS__Vln_plot_muscle-score_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      plot_grid(p1, p2, p3, ncol = 3, align = "h", axis = "tb")
      ggsave(file.path(plot_dir,"7_FPRMS_VlnPlot.pdf"), width=7, height=4.5)
      
      