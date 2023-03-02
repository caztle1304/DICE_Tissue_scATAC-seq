# Final Gene score plots

library(ArchR)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lung CD4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lungCD4 <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')

lungCD4.markerGenes <- scan('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/markerGenes/CD4.txt', what = 'list', sep = '\n')

# All genes
pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/geneScores.pdf')

plotEmbedding(lungCD4, colorBy = 'GeneScoreMatrix',
name = lungCD4.markerGenes,
embedding = 'lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_UMAP',
plotAs = 'points',
size = 0.7)

dev.off()

# Genes we decided to be useful

lungCD4.markerGenes_U <- c('IKZF2', 'CTLA4', 'FOXP3', 'TIGIT', 'ITGA1', 'ZNF683', 'S1PR1', 'GPR183', 'IL7R', 'TCF7', 'SELL', 'PRF1', 'GZMB', 'CRTAM')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD4/geneScores_mostUseful.pdf')

plotEmbedding(lungCD4, colorBy = 'GeneScoreMatrix',
name = lungCD4.markerGenes_U,
embedding = 'lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_UMAP',
plotAs = 'points',
size = 0.7)

dev.off()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lung CD8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lungCD8 <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')

lungCD8.markerGenes <- scan('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/markerGenes/CD8.txt', what = 'list', sep = '\n')

# All genes
pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/geneScores.pdf')

plotEmbedding(lungCD8, colorBy = 'GeneScoreMatrix',
name = lungCD8.markerGenes,
embedding = 'lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_UMAP',
plotAs = 'points',
size = 0.7)

dev.off()

# Genes we decided to be useful

lungCD8.markerGenes_U <- c('ZEB2', 'PRF1', 'KLRG1', 'FGFBP2', 'GNLY', 'FCGR3A', 'TCF7', 'SELL', 'IL7R', 'KLF2', 'CCR7', 'KLRB1', 'JAML', 'ITGAE', 'ITGA1', 'ZNF683')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/lungCD8/geneScores_mostUseful.pdf')

plotEmbedding(lungCD8, colorBy = 'GeneScoreMatrix',
name = lungCD8.markerGenes_U,
embedding = 'lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_UMAP',
plotAs = 'points',
size = 0.7)

dev.off()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Myeloid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myeloid <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k')

myeloid.markerGenes <- scan('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/markerGenes/myeloidMarkerGenes.txt', what = 'list', sep = '\n')

# All genes
pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/geneScores.pdf')

plotEmbedding(myeloid, colorBy = 'GeneScoreMatrix',
name = myeloid.markerGenes,
embedding = 'myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_UMAP',
plotAs = 'points',
size = 0.7)

dev.off()

# Genes we decided to be usfeul

myeloid.markerGenes_U <- c('VCAN', 'S100A8', 'CD14', 'S100A9', 'S100A12', 'IFITM3', 'IFITM1', 'FCGR3A', 'HLA-DQB1', 'HLA-DQA1', 'CCL17',
'C1QC', 'GZMB', 'IRF7', 'WDFY4', 'IRF8', 'LAMP3', 'CCR7', 'CD69', 'HDC')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/myeloid/geneScores_mostUseful.pdf')

plotEmbedding(myeloid, colorBy = 'GeneScoreMatrix',
name = myeloid.markerGenes_U,
embedding = 'myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_UMAP',
plotAs = 'points',
size = 0.7)

dev.off()
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NK <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k')

NK.markerGenes <- scan('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/markerGenes/NKMarkerGenes.txt', what = 'list', sep = '\n')

# All genes
pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/geneScores.pdf')

plotEmbedding(NK, colorBy = 'GeneScoreMatrix',
name = NK.markerGenes,
embedding = 'NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_UMAP',
plotAs = 'points',
size = 0.7)

dev.off()

# Genes we decided to be usfeul

NK.markerGenes_U <- c('CX3CR1', 'TBX21', 'XCL1', 'CXCR3')

pdf('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/finalPlots/NK/geneScores_mostUseful.pdf')

plotEmbedding(NK, colorBy = 'GeneScoreMatrix',
name = NK.markerGenes_U,
embedding = 'NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_UMAP',
plotAs = 'points',
size = 0.7)

dev.off()
