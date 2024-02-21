rm(list = ls())

library(dplyr)
#library(Seurat)
library(ggplot2)
library(patchwork)
#library(SeuratObject)
#library(SingleCellExperiment)
library(readxl)
library(data.table)
#library(paletteer)
#library(RColorBrewer)
#library(scater)
#library(scRNAseq)
#library(SingleR)
#library(SCpubr)


# Set up environment ------------------------------------------------
base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')

output_dir <- file.path(base_dir, 'output/metadata/dev_muscle/distribution')
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

source(file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources/Plot_style_v2.R'))

# color palette
#colors populations
col_celltype <- c('#C70E7BFF', '#FC6882FF',  '#6C6C9DFF',  '#A6E000FF','#1BB6AFFF','#172869FF')
names(col_celltype) <- c('Myogenic Progenitors', 'Myoblast-Myocytes', 'Skeletal mesenchymal', 'Myoblasts', 'Myocytes', 'Postnatal satellite cells')

#colors time points
# col_timepoints <- c('#9E0142FF', '#D53E4FFF', '#F46D43FF', '#F46D43FF', '#FDAE61FF',  '#FEE08BFF',
#                     '#FFFFBFFF',   '#E6F598FF', '#ABDDA4FF', '#66C2A5FF', '#3288BDFF', 
#                     c(rep('#5E4FA2FF', 4)), 'grey')
# names(col_timepoints) <- c("Wk5.0", "Wk6.0" ,  "Wk6.5", "Wk7.0", "Wk7.25",  "Wk7.75" , "Wk9" ,   "Wk12",   "Wk14", "Wk17" ,  "Wk18", 
#                            "Yr7",    "Yr11",   "Yr34" ,  "Yr42", NA)

col_timepoints <- c('#9E0142FF',    '#FDAE61FF',   '#FFFFBFFF',   '#ABDDA4FF', '#66C2A5FF', '#3288BDFF', '#5E4FA2FF')
names(col_timepoints) <- c('Wk5-6', 'Wk6-7', 'Wk7-8' ,'Wk9' ,'Wk12-14', 'Wk17-18', 'Postnatal')


# Load metadata ------------------------------------------------
P3F1 <- readRDS(file.path(base_dir, "output/metadata/dev_muscle/metadata/P3F1 FP-RMS metadata.rds"))
P7F1 <- readRDS(file.path(base_dir, "output/metadata/dev_muscle/metadata/P7F1 FP-RMS metadata.rds"))
FNRMS <- readRDS(file.path(base_dir, "output/metadata/dev_muscle/metadata/FN-RMS metadata.rds"))

metadata <- list(P3F1, P7F1, FNRMS)

# subset to columns of interest
metadata_red <- list()
for (i in seq_along(metadata)){
  metadata_red[[i]] <- metadata[[i]][ ,c('name', 'SingleR.cell.labels', 'SingleR.labels')]
}

# concatenate
metadata <- rbind(metadata_red[[1]], metadata_red[[2]], metadata_red[[3]])
desired_order <- c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3',  #P3F1
                   'aRMS-4',  'aRMS-5', 'KFR', 'Mast118', #P3F1
                   'MSK82489', 'Rh4',  'Rh41',   'RMS', 
                   'SJRHB013759_A1',  'SJRHB013759_A2', 
                   'SJRHB013759_X14','SJRHB013759_X15', 
                   'Mast95', 'MSK72117', 'MSK72117_SC', #P7F1
                   'SJRHB010468_D1',  'SJRHB010468_X1', 'SJRHB013757_D2', 'SJRHB013757_X1', #P7F1
                   'SJRHB031320_D1',   'SJRHB031320_X1', 'SJRHB046156_A1', 'SJRHB046156_X1', #P7F1
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
                   'SJRHB030680_X1', 'SJRHB049189_D1', 'SJRHB049189_X1')

# Plot distribution of cell types  ------------------------------------------------

# Calculate frequency 
metadata_freq <- metadata %>% 
  dplyr::group_by(name, SingleR.cell.labels) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n)) 
metadata_freq

# Plot
metadata_freq$name <- factor(metadata_freq$name, levels = desired_order)

gg_barplot_style(metadata_freq, col_celltype) +
  geom_bar(mapping = aes(x = name, y = freq, fill = SingleR.cell.labels), stat = "identity", width = 0.9, color="black") +
  guides(fill = guide_legend(ncol = 1, title = 'Cell Type')) +
  labs(x = 'Sample', y = 'Percentage') 
ggsave(file.path(output_dir, "4_Barplot_composition_RMS.pdf"), width=18, height=4, dpi=300)


# Plot distribution of time points ------------------------------------------------

# Calculate frequency 
metadata_freq <- metadata %>% 
  dplyr::group_by(name, SingleR.labels) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n)) 
metadata_freq

# Plot
metadata_freq$name <- factor(metadata_freq$name, levels = desired_order)

gg_barplot_style(metadata_freq, col_timepoints) +
  geom_bar(mapping = aes(x = name, y = freq, fill = SingleR.labels), stat = "identity", width = 0.9, color="black") +
  guides(fill = guide_legend(ncol = 1, title = 'Cell Type')) +
  labs(x = 'Sample', y = 'Percentage') 
ggsave(file.path(output_dir, "4_Barplot_composition_RMS_timepoints.pdf"), width=18, height=4, dpi=300)



