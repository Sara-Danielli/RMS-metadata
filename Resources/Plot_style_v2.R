library(dplyr)
library(Seurat)
library(ggplot2)

## color palette

# color cluster names FPRMS/FNRMS
col_cluster_names_aggregate <- c("#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#8E0152FF",'#66C5CCFF' , '#D8D8E0FF', '#B497E7FF')
names(col_cluster_names_aggregate) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis', 'IFN')

# color cluster names RMS atlas
col_cluster_names_aggregate_integrated <- c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF", "#8E0152FF", '#D8D8E0FF')
names(col_cluster_names_aggregate_integrated) <- c('Progenitor','TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated', 'Differentiated', 'Apoptosis')

# Color model
col_model <- c("#009B9EFF","#A7D3D4FF",  "#E4C1D9FF","#C75DABFF")
names(col_model) <- c("Patient", "O-PDX", "Primary culture", "Cell line")

# Color subtype
col_subtype <- c('#D3A2C2FF', '#95CECFFF')

# Color origin
col_origin <- paletteer::paletteer_d("ggthemes::excel_Slice", n=4)
names(col_origin) <- c("Wei et al.", "Patel et al.", "Danielli et al.", "Weng et al.")

# Color name
col_aRMS <- paletteer::paletteer_c("ggthemes::Blue-Teal", n = 27)
col_eRMS <- paletteer::paletteer_c("ggthemes::Purple", n = 47)
col_name <- c(col_aRMS, col_eRMS)
names(col_name) <- c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3', 
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
                    'SJRHB030680_X1', 'SJRHB049189_D1', 'SJRHB049189_X1')

# Color subtype
col_subtype <- c('#D3A2C2FF', '#95CECFFF')



## Style plots

# Barplot function
theme_ggplot = theme(legend.position = "none",
                     plot.title = element_text(hjust=0.5, face="bold"),
                     panel.border = element_blank(),
                     plot.background = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1, colour="black"),
                     axis.text.y = element_text(size=12, colour="black"),
                     axis.title=element_text(size=12),
                     strip.text = element_text(size = 13, face = "bold"))

theme_ggplot_legend = theme(plot.title = element_text(hjust=0.5, face="bold"),
                            panel.border = element_blank(),
                            plot.background = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1, colour="black"),
                            axis.text.y = element_text(size=12, colour="black"),
                            axis.title=element_text(size=12),
                            strip.text = element_text(size = 13, face = "bold"))

## Style plots

theme_vln <- theme(panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                   axis.text.y = element_text(size=14, colour="black"),
                   axis.title=element_text(size=14),
                   plot.title = element_text(size=14, face="bold")) 


theme_cellstate_plot <-  theme(panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               plot.background = element_blank(),
                               panel.background = element_blank(),
                               axis.line = element_line(colour = "black"),
                               axis.title=element_text(size=16),
                               axis.text.x = element_text(size=16, vjust = 0.5, colour="black"),
                               axis.text.y = element_text(size=16, colour="black"),
                               plot.title = element_text(hjust = 0.5, angle = 0, size = 16, face = "bold", vjust = 1),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 15, face = "bold"),
                               legend.box.background = element_rect(color="black", size=1),
                               strip.text = element_text(size = 15, face = "bold"))

# Barplot function
plot_bar <- function(seurat_obj, x_var, y_var, colors){
  ggplot(seurat_obj@meta.data, aes(x_var, fill = y_var)) +
    scale_fill_manual(values = colors) + 
    geom_bar(position = "fill", color="black") +
    labs (y='Proportion', x='') + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=14, colour="black"),
          axis.title=element_text(size=14),
          legend.text = element_text(size = 14),
          legend.title = element_blank()) 
}
