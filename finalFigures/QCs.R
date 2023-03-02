############    -  Plotting final QC metrics  -    ############

# ---
# Author: Angel Adrian De la Cruz Castillo
# Date: 2023-03-02
# ---

### -------------------------- Description -------------------------- ###
# This script creates all final QC metrics used in paper: TSS vs Fragments plot, fragment size distribution, TSS enrichment and a comined 
# violin/boxplot containing prediction scores for each cell type resulting from scRNA-seq integration


library(ArchR)

plotTSSvsFrags <- function(project, xintercept, yintercept) {
df <- getCellColData(project, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1],
    y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(2.5, max(df[,1] * 1.05)),
    ylim = c(0, max(df[,2] * 1.05))
) + geom_hline(yintercept = yintercept, lty = "dashed") + geom_vline(xintercept = log10(xintercept), lty = "dashed")

return(p)
}

plotPredictionScorePerCelltype <- function(archrproj, suffix)
{
  ## Plots violin plot with boxplot for scoring of each cell type
cellColData <- getCellColData(archrproj, select = paste0(c('predictedGroup_', 'predictedScore_'), suffix))
cellColData <- as.data.frame(cellColData)

p <- ggplot(cellColData, aes(x = get(paste0('predictedGroup_', suffix)), y = get(paste0('predictedScore_', suffix)), fill = get(paste0('predictedGroup_', suffix)))) +
     geom_violin(width = 0.5, alpha = 0.2) +
     geom_boxplot(width = 0.1) +
     theme_bw() +
     ggtitle('Prediction scores by group') +
     theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
     xlab('Cell type') + ylab('Prediction score')
return(p)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lung CD4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lungCD4 <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')
lungCD4.unfiltered <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/unfiltered_projects/lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')

# Fragment size distribution
pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/fragmentSizeDistribution.pdf')

plotFragmentSizes(lungCD4, groupBy = 'Sample')

dev.off()

# TSS enrichment

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/TSSEnrichment.pdf')

plotTSSEnrichment(lungCD4, groupBy = 'Sample')

dev.off()

# Violin plot with integration scores per cell type

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/integrationScores_perCellType.pdf')

plotPredictionScorePerCelltype(lungCD4, 'Un')

dev.off()

# TSSvsFrags

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/TSSvsFrags.pdf')

plotTSSvsFrags(lungCD4.unfiltered, 3162, 10)

dev.off()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lung CD8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lungCD8 <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')
lungCD8.unfiltered <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/unfiltered_projects/lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/fragmentSizeDistribution.pdf')

plotFragmentSizes(lungCD8, groupBy = 'Sample')

dev.off()

# TSS enrichment

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/TSSEnrichment.pdf')

plotTSSEnrichment(lungCD8, groupBy = 'Sample')

dev.off()

# Violin plot with integration scores per cell type

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/integrationScores_perCellType.pdf')

plotPredictionScorePerCelltype(lungCD8, 'Un')

dev.off()

# TSSvsFrags

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/TSSvsFrags.pdf')

plotTSSvsFrags(lungCD8.unfiltered, 3162, 10)

dev.off()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NK <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k')
NK.unfiltered <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/unfiltered_projects/NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/fragmentSizeDistribution.pdf')

plotFragmentSizes(NK, groupBy = 'Sample')

dev.off()

# TSS enrichment

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/TSSEnrichment.pdf')

plotTSSEnrichment(NK, groupBy = 'Sample')

dev.off()

# Violin plot with integration scores per cell type

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/integrationScores_perCellType.pdf')

plotPredictionScorePerCelltype(NK, 'Un')

dev.off()

# TSSvsFrags

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/TSSvsFrags.pdf')

plotTSSvsFrags(NK.unfiltered, 3162, 8)

dev.off()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Myeloid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###QCs will be calculated for ALL cells passing QCs, meaning MAST cell will be included here
myeloid <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k')
myeloid.unfiltered <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/unfiltered_projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/fragmentSizeDistribution.pdf')

plotFragmentSizes(myeloid, groupBy = 'Sample')

dev.off()

# TSS enrichment

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/TSSEnrichment.pdf')

plotTSSEnrichment(myeloid, groupBy = 'Sample')

dev.off()

# Violin plot with integration scores per cell type

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/integrationScores_perCellType.pdf')

plotPredictionScorePerCelltype(myeloid, 'Un')

dev.off()

# TSSvsFrags

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/TSSvsFrags.pdf')

plotTSSvsFrags(myeloid.unfiltered, 3162, 8)

dev.off()
