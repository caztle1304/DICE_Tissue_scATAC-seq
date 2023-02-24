library(ArchR)

myeloid <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k')

cells2keep <- myeloid$cellNames[which(myeloid$myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_1.7 != 'C16')]

myeloid.subset <- subsetArchRProject(ArchRProj = myeloid, cells = cells2keep, outputDirectory = '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_noMAST')
