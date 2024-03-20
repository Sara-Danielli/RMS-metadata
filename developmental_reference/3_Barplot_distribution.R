rm(list = ls())

library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(data.table)
library(ggsankey)


# Set up environment ------------------------------------------------
base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')

output_dir <- file.path(base_dir, 'output/metadata/dev_muscle/distribution')
if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}

source(file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources/Plot_style_v2.R'))

# color palette
#colors populations
#col_celltype_pal2 <- c('#C3A016FF', '#C3D878FF', '#58A787FF', '#8EBACDFF', '#246893FF', '#163274FF', '#0C1F4BFF')
#col_celltype_pal0 <- c('#C70E7BFF', '#FC6882FF',  '#A6E000FF','#1BB6AFFF', '#6C6C9DFF', '#172869FF')
#col_celltype_pal0b <- c('#C70E7BFF', '#FC6882FF', '#1BB6AFFF', '#6C6C9DFF', '#172869FF', 'black')
#col_celltype_pal1 <- c('#62205FFF', '#D0A8C8FF', '#DCEAEAFF', '#7CC3C9FF', '#172869FF', '#FBF5A3FF')
col_celltype <- c('#62205FFF', '#D0A8C8FF', '#DCEAEAFF', '#7CC3C9FF', '#172869FF', '#C3D878FF')
#col_celltype_pal3 <- c('#CF5483FF', '#FE8BB2FF', '#B6CAC1FF', '#98E2E3FF' ,'#57ADBAFF', '#5F4E3CFF')
#col_celltype_pal4 <- paletteer::paletteer_d("MoMAColors::Rattner")[1:6]
names(col_celltype) <- c('Skeletal mesenchymal', 'Myogenic Progenitors', 'Myoblast-Myocytes', 'Myoblasts', 'Myocytes', 'Postnatal satellite cells')

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

metadata <- list(FNRMS, P7F1, P3F1)


# subset to columns of interest
metadata_red <- list()
for (i in seq_along(metadata)){
  metadata_red[[i]] <- metadata[[i]][ ,c('name', 'SingleR.cell.labels', 'SingleR.labels')]
}

# concatenate
metadata <- rbind(metadata_red[[1]], metadata_red[[2]], metadata_red[[3]])
desired_order <- c('20696','21202', '29806', 
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
                   'SJRHB030680_X1', 'SJRHB049189_D1', 'SJRHB049189_X1',
                   'Mast95', 'MSK72117', 'MSK72117_SC', #P7F1
                   'SJRHB010468_D1',  'SJRHB010468_X1', 'SJRHB013757_D2', 'SJRHB013757_X1', #P7F1
                   'SJRHB031320_D1',   'SJRHB031320_X1', 'SJRHB046156_A1', 'SJRHB046156_X1', #P7F1
                   '20082',  'aRMS-1',  'aRMS-2', 'aRMS-3',  #P3F1
                   'aRMS-4',  'aRMS-5', 'KFR', 'Mast118', #P3F1
                   'MSK82489', 'Rh4',  'Rh41',   'RMS', 
                   'SJRHB013759_A1',  'SJRHB013759_A2', 
                   'SJRHB013759_X14','SJRHB013759_X15'
)


# Plot distribution of cell types  ------------------------------------------------

# Calculate frequency 
metadata_freq <- metadata %>% 
  dplyr::group_by(name, SingleR.cell.labels) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n)) 
metadata_freq

# Plot
metadata_freq$name <- factor(metadata_freq$name, levels = desired_order)
metadata_freq$SingleR.cell.labels <- factor(metadata_freq$SingleR.cell.labels,
                                            levels = c('Skeletal mesenchymal', 'Myogenic Progenitors', 'Myoblast-Myocytes', 'Myoblasts', 'Myocytes', 'Postnatal satellite cells'))


gg_barplot_style(metadata_freq, col_celltype) +
  geom_bar(mapping = aes(x = name, y = freq, fill = SingleR.cell.labels), stat = "identity", width = 0.9, color="black") +
  guides(fill = guide_legend(ncol = 1, title = 'Cell Type')) +
  labs(x = 'Sample', y = 'Percentage') 
ggsave(file.path(output_dir, "1_Barplot_composition_RMS_pal0b.pdf"), width=18, height=4, dpi=300)


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
ggsave(file.path(output_dir, "2_Barplot_composition_RMS_timepoints.pdf"), width=18, height=4, dpi=300)



# Load .rds objects to plot alluvial plots ------------------------------------------------
P7F1_sc <- readRDS(file.path(base_dir, "write/FPRMS_PAX7FOXO1_final_20240130.rds"))
P3F1_sc <- readRDS(file.path(base_dir, "write/FPRMS_PAX3FOXO1_final_20240130.rds"))
FNRMS_sc <- readRDS(file.path(base_dir, "write/FNRMS_final_20240130.rds"))

# add SingleR label
P7F1_sc <- AddMetaData(P7F1_sc, P7F1, col.name = 'SingleR.cell.labels')
P3F1_sc <- AddMetaData(P3F1_sc, P3F1, col.name = 'SingleR.cell.labels')
FNRMS_sc <- AddMetaData(FNRMS_sc, FNRMS, col.name = 'SingleR.cell.labels')

# extract metadata for Sankey plots 
P7F1_metadata <- P7F1_sc@meta.data[ , c('Cluster assignment', 'SingleR.cell.labels')]
P3F1_metadata <- P3F1_sc@meta.data[ , c('Cluster assignment', 'SingleR.cell.labels')]
FNRMS_metadata <- FNRMS_sc@meta.data[ , c('Cluster assignment', 'SingleR.cell.labels')]


# create Sankey plots ----------------------------------------------------------
metadata <- list(P7F1_metadata, P3F1_metadata, FNRMS_metadata)
names(metadata) <- c('PAX7FOXO1', 'PAX3FOXO1', 'FNRMS')

df <- list()

for (i in seq_along(metadata)){
  # transform data for plot
  df[[i]] <- metadata[[i]] %>%
    make_long(`Cluster assignment`, SingleR.cell.labels)

  # Reorder factors
  desired_order_node <- rev(c('Skeletal mesenchymal', 'Myogenic Progenitors', 'Myoblasts', 
                              'Myoblast-Myocytes',  'Myocytes', 'Postnatal satellite cells',
                              'Progenitor', 'Proliferative', 'Ground', 'Differentiated', 
                              'Neuronal', 'Apoptosis', 'IFN'))
  
  df[[i]]$node <- factor(df[[i]]$node,levels = desired_order_node)
  df[[i]]$next_node <- factor(df[[i]]$next_node,levels = desired_order_node)
}
  
  # Chart 
for (i in seq_along(metadata)){
  ggplot(df[[i]], aes(x = x
                        , next_x = next_x
                        , node = node
                        , next_node = next_node
                        , fill = factor(node)
                        , label = node)
  ) +geom_sankey(flow.alpha = 0.3
                 , node.color = "black"
                 ,show.legend = FALSE) +
    geom_sankey_label(size = 3.5, color = 1, fill = "white", 
                      hjust = -0.1) +
    theme_bw() + 
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) + 
    scale_fill_manual(values = c(col_cluster_names_aggregate, col_celltype)) + 
    labs(title = "FN-RMS")
  ggsave(file.path(output_dir, paste0("3_SankeyPlot_", names(metadata)[i], ".pdf")), width=5, height=5)
 
}


# Plot UMAP plots with SingleR labels -----------------------------------------
p1 <- DimPlot(FNRMS_sc, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col_celltype, 
              label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + ggtitle('FN-RMS')
p2 <- DimPlot(P7F1_sc, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col_celltype, 
              label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + NoLegend() + ggtitle('FP-RMS (PAX7::FOXO1)')
p3 <- DimPlot(P3F1_sc, reduction = "umap_rpca", group.by = 'SingleR.cell.labels', cols = col_celltype, 
              label=F, pt.size = 3, raster=TRUE, raster.dpi = c(1012, 1012)) + ggtitle('FP-RMS (PAX3::FOXO1')
p1| p2 | p3
ggsave(file.path(output_dir, "4_UMAP_cell_projections.pdf"), width=12.5, height=4, dpi=300)


### SkM. mesen population split by primary/recurrence
Idents(FNRMS_sc) <- 'model'
FNRMS_sc_subset <- subset(FNRMS_sc, idents = "Patient")
Idents(FNRMS_sc_subset) <- 'treatment'
FNRMS_sc_subset <- subset(FNRMS_sc_subset, idents = "NA", invert = T)

FNRMS_sc_subset$SingleR.cell.labels <- factor(FNRMS_sc_subset$SingleR.cell.labels,
                                            levels = c('Skeletal mesenchymal', 'Myogenic Progenitors', 'Myoblast-Myocytes', 'Myoblasts', 'Myocytes', 'Postnatal satellite cells'))

### SkM. mesen population split by primary/recurrence
FNRMS_sc_meta <- FNRMS_sc_subset@meta.data %>% 
  group_by(name, PatientID, treatment, SingleR.cell.labels) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n)) 

ggplot(FNRMS_sc_meta, aes(x=treatment, y=freq)) +
  geom_violin(aes(fill = SingleR.cell.labels)) + 
  geom_boxplot(width = 0.07) +
  geom_point(fill = 'black', aes(group = PatientID),size=2,shape=21, position = position_dodge(0.2)) +
  geom_line(aes(group = PatientID), position = position_dodge(0.2), linetype = "dashed") +
  theme_vln_plot +
  labs (y='Proportion', x='') + 
  scale_fill_manual(values=col_celltype) + 
  facet_wrap(. ~ SingleR.cell.labels, nrow = 1, scales = "free") +
  stat_compare_means(aes(label = after_stat(p.format)), 
                     method = 't.test', 
                     #ref.group = 'FN-RMS',
                     size = 4, 
                     label.y.npc = 0.91, 
                     label.x.npc = 0.4)  + NoLegend() +
  scale_y_continuous(labels = scales::percent) 
ggsave(file.path(output_dir, "6_VlnPlot_FNRMS.pdf"), width=15, height=3.5, dpi=300)



### Tumor population split by primary/recurrence
FNRMS_sc_meta <- FNRMS_sc_subset@meta.data %>% 
  group_by(name, PatientID, treatment, `Cluster assignment`) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n)) 

ggplot(FNRMS_sc_meta, aes(x=treatment, y=freq)) +
  geom_violin(aes(fill = `Cluster assignment`)) + 
  geom_boxplot(width = 0.07) +
  geom_point(fill = 'black', aes(group = PatientID),size=2,shape=21, position = position_dodge(0.2)) +
  geom_line(aes(group = PatientID), position = position_dodge(0.2), linetype = "dashed") +
  theme_vln_plot +
  labs (y='Proportion', x='') + 
  scale_fill_manual(values=col_cluster_names_aggregate) + 
  facet_wrap(. ~ `Cluster assignment`, nrow = 1, scales = "free") +
  stat_compare_means(aes(label = after_stat(p.format)), 
                     method = 't.test', 
                     #ref.group = 'FN-RMS',
                     size = 4, 
                     label.y.npc = 0.91, 
                     label.x.npc = 0.4)  + NoLegend() +
  scale_y_continuous(labels = scales::percent) 
ggsave(file.path(output_dir, "6_VlnPlot_FNRMS_tumor_pop.pdf"), width=10, height=3.5, dpi=300)


