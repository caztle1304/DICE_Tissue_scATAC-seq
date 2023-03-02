############    -  UMAP plotting  -    ############

# ---
# Author: Angel Adrian De la Cruz Castillo
# Date: 2023-03-02
# ---

### -------------------------- Description -------------------------- ###
# This script was used to plot UMAPs colored by different groups: unbiased clustering, results of integration with scRNA-seq and prediction 
# scores of integration with scRNA-seq

library(ArchR)
library(stringr)

plotUMAP <- function(df, color_by, colors, size, title){
  p <- ggplot(data = df, aes(x = .data[["Dimension_1"]], y = .data[["Dimension_2"]], color=.data[[color_by]])) + geom_point(size = size) + scale_color_manual(values = colors) + theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) + ggtitle(title) + guides(colour = guide_legend(override.aes = list(size = 2)))
  return(p)
}

plotUMAP_empty <- function(df, color_by, colors, size){
  p <- ggplot(data = df, aes(x = .data[["Dimension_1"]], y = .data[["Dimension_2"]], color=.data[[color_by]])) + geom_point(size = size) + scale_color_manual(values = colors) + theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'none',
    panel.border = element_blank()

  )
  return(p)
}


retrieveUMAP <- function(project, embedding, dimReduction){
metrics <- getCellColData(project)
emb <- getEmbedding(project, embedding = embedding)
cells <- intersect(rownames(metrics),rownames(emb))
metrics <- metrics[order(cells),]
emb <- emb[order(cells),]
metrics["Dimension_1"] <- emb[,paste0(dimReduction,"#UMAP_Dimension_1")]
metrics["Dimension_2"] <- emb[,paste0(dimReduction,"#UMAP_Dimension_2")]
metrics <- as.data.frame(metrics)
return(metrics)

}

plotUMAP_continuous <- function(df, color_by, size, colors, title){
  p <- ggplot(data = df, aes(x = .data[["Dimension_1"]], y = .data[["Dimension_2"]], color=.data[[color_by]])) + geom_point(size = size) +
  scale_colour_gradientn(colours = colors) + theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank()
  ) + ggtitle(title)
  return(p)
}

plotUMAP_continuous_empty <- function(df, color_by, colors, size){
  p <- ggplot(data = df, aes(x = .data[["Dimension_1"]], y = .data[["Dimension_2"]], color=.data[[color_by]])) + geom_point(size = size) +
  scale_colour_gradientn(colours = colors) + theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'none',
    panel.border = element_blank()

  )
  return(p)
}


stallion = c("C1"="#D51F26","C2"="#272E6A","C3"="#208A42","C4"="#89288F","C5"="#F47D2B", "C6"="#FEE500","C7"="#8A9FD1","C8"="#C06CAB","C9"="#D8A767",
               "C10"="#90D5E4", "C11"="#89C75F","C12"="#F37B7D","C13"="#9983BD","C14"="#D24B27","C15"="#3BBCA8", "C16"="#6E4B9E","C17"="#0C727C", "C18"="#7E1416","C19"="#E6C2DC","C20"="#3D3D3D")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lung CD4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lungCD4 <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')

lungCD4.coords <- retrieveUMAP(lungCD4, 'lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_UMAP', 'lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony')

## Naming properly manual annotations
lungCD4.coords$annotations_12_23_22 <- str_replace(lungCD4.coords$annotations_12_23_22, 'Intermediate_population', 'TCM_TRM')

## Naming properly integration results
lungCD4.coords$predictedGroup_Un_Round2 <- str_replace(lungCD4.coords$predictedGroup_Un_Round2, 'Bonafide_CD4_TRMs', 'TRM')
lungCD4.coords$predictedGroup_Un_Round2 <- str_replace(lungCD4.coords$predictedGroup_Un_Round2, 'TRM_w_circulatingFeatures', 'TCM_TRM')


colors <- c('TCM_TRM' = '#00BFFF', 'CD4_CTL' = '#E60000', 'TCM' = '#BFEFFF', 'TRM' = '#3333FF', 'TREG' = '#319272') # For manual annotations

colors2 <- c('TCM_TRM' = '#00BFFF', 'Conventional_CD4_CTL' = '#E60000', 'TCM' = '#BFEFFF', 'TRM' = '#3333FF', 'TREG' = '#319272', 'Non_conventional_CD4_CTL' = '#FF9933') # For integration results

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/unbiased_clustering.pdf')

plotUMAP(lungCD4.coords, 'lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.8',
stallion[seq(1, length(unique(lungCD4$lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.8)))], 0.7, 'Unbiased clustering, resolution 0.8')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/unbiased_clustering_noLabels.pdf')

plotUMAP_empty(lungCD4.coords, 'lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.8',
stallion[seq(1, length(unique(lungCD4$lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.8)))], 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/manual_annotations_December_2022.pdf')

plotUMAP(lungCD4.coords, 'annotations_12_23_22', colors, 0.7, 'Manual Annotations December 23 2022')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/manual_annotations_December_2022_noLabels.pdf')

plotUMAP_empty(lungCD4.coords, 'annotations_12_23_22', colors, 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/integrationResults.pdf')

plotUMAP(lungCD4.coords, 'predictedGroup_Un_Round2', colors2, 0.7, 'Label transferring with scRNA-seq integration')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/integrationResults_noLabels.pdf')

plotUMAP_empty(lungCD4.coords, 'predictedGroup_Un_Round2', colors2, 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/predictionScores_integration.pdf')

plotUMAP_continuous(lungCD4.coords, 'predictedScore_Un_Round2', 0.7, ArchRPalettes$solarExtra, 'Prediction scores from integration')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/predictionScores_integration_noLabels.pdf')

plotUMAP_continuous_empty(lungCD4.coords, 'predictedScore_Un_Round2', ArchRPalettes$solarExtra, 0.7)

dev.off()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lung CD8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lungCD8 <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')

lungCD8.coords <- retrieveUMAP(lungCD8, 'lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_UMAP', 'lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony')


colors <- c('GZMK_expressing' = '#4DA6FF', 'Effector_memory' = '#CC99FF', 'TRM' = '#8000FF', 'TCM' = '#FF1AB3', 'MAIT' = '#676765') # For manual annotations

colors2 <- c('GZMK_expressing' = '#4DA6FF', 'Effector_memory' = '#CC99FF', 'TRM' = '#8000FF', 'TCM' = '#FF1AB3', 'MAIT' = '#676765', 'Inntate_receptor_expressing' = '#39AC39') # For integration results

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/unbiased_clustering.pdf')

plotUMAP(lungCD8.coords, 'lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.4',
stallion[seq(1, length(unique(lungCD8$lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.4)))], 0.7, 'Unbiased clustering, resolution 0.4')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/unbiased_clustering_noLabels.pdf')

plotUMAP_empty(lungCD8.coords, 'lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.4',
stallion[seq(1, length(unique(lungCD8$lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.4)))], 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/manual_annotations_December23_2022.pdf')

plotUMAP(lungCD8.coords, 'annotations_12_23_22', colors, 0.7, 'Manual Annotations December 23 2022')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/manual_annotations_December23_2022_noLabels.pdf')

plotUMAP_empty(lungCD8.coords, 'annotations_12_23_22', colors, 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/integrationResults.pdf')

plotUMAP(lungCD8.coords, 'predictedGroup_Un_Round2', colors2, 0.7, 'Label transferring with scRNA-seq integraiton')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/integrationResults_noLabels.pdf')

plotUMAP_empty(lungCD8.coords, 'predictedGroup_Un_Round2', colors2, 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/predictionScores_integration.pdf')

plotUMAP_continuous(lungCD8.coords, 'predictedScore_Un_Round2', 0.7, ArchRPalettes$solarExtra, 'Prediction scores from integration')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/predictionScores_integration_noLabels.pdf')

plotUMAP_continuous_empty(lungCD8.coords, 'predictedScore_Un_Round2', ArchRPalettes$solarExtra, 0.7)

dev.off()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NK <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k')

NK.coords <- retrieveUMAP(NK, 'NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_UMAP', 'NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony')


# Naming properly manual annotations
NK.coords$annotations_12_23_22 <- str_replace(NK.coords$annotations_12_23_22, 'CD16Positive', 'CD56dim/CD16+')
NK.coords$annotations_12_23_22 <- str_replace(NK.coords$annotations_12_23_22, 'CD16Negative', 'CD56bright/CD16-')

# Naming properly integration results
NK.coords$predictedGroup_Un_Round2 <- str_replace(NK.coords$predictedGroup_Un_Round2, 'CD16Positive', 'CD56dim/CD16+')
NK.coords$predictedGroup_Un_Round2 <- str_replace(NK.coords$predictedGroup_Un_Round2, 'CD16Negative', 'CD56bright/CD16-')

colors <- c('CD56dim/CD16+' = '#FFD700', 'CD56bright/CD16-' = '#E87600')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/unbiased_clustering.pdf')

plotUMAP(NK.coords, 'NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_0.4', stallion[seq(length(unique(NK$NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_0.4)))], 0.7, 'Unbiased clustering, resolution 0.4')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/unbiased_clustering_noLabels.pdf')

plotUMAP_empty(NK.coords, 'NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_0.4', stallion[seq(length(unique(NK$NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_0.4)))], 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/manual_annotations_December23_2022.pdf')

plotUMAP(NK.coords, 'annotations_12_23_22', colors, 0.7, 'Manual Annotations December 23 2022')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/manual_annotations_December23_2022_noLabels.pdf')

plotUMAP_empty(NK.coords, 'annotations_12_23_22', colors, 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/integrationResults.pdf')

plotUMAP(NK.coords, 'predictedGroup_Un_Round2', colors, 0.7, 'Label transferring with scRNA-seq integration')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/integrationResults_noLabels.pdf')

plotUMAP_empty(NK.coords, 'predictedGroup_Un_Round2', colors, 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/predictionScores_integration.pdf')

plotUMAP_continuous(NK.coords, 'predictedScore_Un_Round2', 0.7, ArchRPalettes$solarExtra, 'Prediction scores from integration')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/predictionScores_integration_noLabels.pdf')

plotUMAP_continuous_empty(NK.coords, 'predictedScore_Un_Round2', ArchRPalettes$solarExtra, 0.7)

dev.off()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Myeloid No MAST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myeloid <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_noMAST')

myeloid.coords <- retrieveUMAP(myeloid, 'myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_UMAP', 'myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony')

colors <- c('Classical_monocytes' = '#B22222', 'Tissue_resident_Macrophages' = '#0052CC', 'pDCs' = '#BFBFBF', 'mregDC' = '#9EC026', 'Non_classical_Monocytes' = '#FF8C00',
          'Macrophage_cDC2_mixture' = '#80BFFF', 'cDC2s' = '#EFA9A9', 'cDC1s' = '#803E96')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/integrationResults.pdf')

plotUMAP(myeloid.coords, 'predictedGroup_Un', colors, 0.7, 'Label transferring with scRNA-seq integration')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/integrationResults_noLabels.pdf')

plotUMAP_empty(myeloid.coords, 'predictedGroup_Un', colors, 0.7)

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/predictionScores_integration.pdf')

plotUMAP_continuous(myeloid.coords, 'predictedScore_Un', 0.7, ArchRPalettes$solarExtra, 'Prediction scores from integration')

dev.off()

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/predictionScores_integration_noLabels.pdf')

plotUMAP_continuous_empty(myeloid.coords, 'predictedScore_Un', ArchRPalettes$solarExtra, 0.7)

dev.off()
