############    -  Removing MAST cells cluster -    ############

# ---
# Author: Angel Adrian De la Cruz Castillo
# Date: 2023-02-15
# ---

### -------------------------- Description -------------------------- ###
# This script will remove cluster of MAST cells from myeloid cells and create a new ArchR project to avoid overwriting

library(ArchR)

myeloid <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k')

cells2keep <- myeloid$cellNames[which(myeloid$myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_1.7 != 'C16')]

myeloid.subset <- subsetArchRProject(ArchRProj = myeloid, cells = cells2keep, outputDirectory = '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_noMAST')
