library(ArchR)
set.seed(1)
addArchRGenome('hg38')

createArchRProject <- function(inputFiles, sampleNames, outputDirectory, frags, TSS)
{
  ArrowFiles <- createArrowFiles(inputFiles = inputFiles, sampleNames = sampleNames,
                  filterTSS = TSS, filterFrags = frags,
                addTileMat = TRUE, addGeneScoreMat = TRUE,
              subThreading = FALSE)
  doubScores <- addDoubletScores(input = ArrowFiles, k = 10, knnMethod = 'UMAP',
                LSIMethod = 1, threads = 1)

  Project <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = outputDirectory, copyArrows = TRUE)

  Project <- filterDoublets(Project)

  saveArchRProject(Project, load = TRUE)
  return(Project)
}

## Lung CD4
lungCD4.sample <- 'R24P_Hu_2_CD4_4D_ATAC'
lungCD4.inputFile <- '/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/CD4_CD8_cellrangerV2/manualNV047/count/R24P_Hu_2_CD4_4D_ATAC/outs/fragments.tsv.gz'
lungCD4.outputDirectory <- '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/unfiltered_projects/lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k'
lungCD4.project <- createArchRProject(lungCD4.inputFile,lungCD4.sample, lungCD4.outputDirectory, 0, 0)

## LungCD8 
lungCD8.sample <- 'R24P_Hu_3_CD8_4D_ATAC'
lungCD8.inputFile <- '/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/CD4_CD8_cellrangerV2/manualNV047/count/R24P_Hu_3_CD8_4D_ATAC/outs/fragments.tsv.gz'
lungCD8.outputDirectory <- '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/unfiltered_projects/lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k'
lungCD8.project <- createArchRProject(lungCD8.inputFile, lungCD8.sample, lungCD8.outputDirectory, 0, 0)

## Myeloid 
myeloid2.sample <- 'R24DICE_Hu_M2_mye_ATAC' 
myeloid2.inputFile <- '/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/NV068/NV068_cellRangerV2/NV068/count/R24DICE_ATAC/R24DICE_Hu_M2_mye_ATAC/outs/fragments.tsv.gz'
myeloid2.outputDirectory <- '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/unfiltered_projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k'
myeloid2.project <- createArchRProject(myeloid2.inputFile, myeloid2.sample, myeloid2.outputDirectory, 0, 0)

## NK
NK1.sample <- 'R24DICE_Hu_N1_NK_ATAC'
NK1.inputFile <- '/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/NV068/NV068_cellRangerV2/NV068/count/R24DICE_ATAC/R24DICE_Hu_N1_NK_ATAC/outs/fragments.tsv.gz'
NK1.outputDirectory <- '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/unfiltered_projects/NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k'
NK1.project <- createArchRProject(NK1.inputFile, NK1.sample, NK1.outputDirectory, 0, 0)
