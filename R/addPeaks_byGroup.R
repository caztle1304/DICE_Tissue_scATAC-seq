library(ArchR)
library(stringr)
library(dplyr)
library(argparse)

args <- ArgumentParser()

args$add_argument('--Project', default = NULL, type = 'character', help = 'Path to ArchR project')
args$add_argument('--groupColName', default = NULL, type = 'character', help = 'Name of column to use for grouping')
args$add_argument('--outDir', default = NULL, type = 'character', help = 'Path to output directory where peaks will be written')
args$add_argument('--cutOff', default = 0.1, type = 'double', help = 'q-value cutoff for peak calling')
args$add_argument('--chromSizes', default = NULL, type = 'character', help = 'Path to chromSizes file for bedToBigBed')
args$add_argument('--bedSort', default = NULL, type = 'character', help = 'Path to bedSort executable')
args$add_argument('--bedToBigBed', default = NULL, type = 'character', help = 'Path to bedToBigBed executable')

args <- args$parse_args()

getPeaks <- function(project)
{
  peaks <- getPeakSet(project)
  peakGroups <- names(peaks)
  names(peaks) <- NULL
  peaks <- as.data.frame(peaks)
  peaks$group <- peakGroups
  return(peaks)
}

writePeaksAndRunBigBed <- function(outDir, peaks, name, bedSort, bedToBigBed, chromSizes)
{
  peaks <- as.data.frame(peaks)
  write.csv(peaks, file = paste0(outDir, 'peaksGroupedBy_', name, '_full.csv'), row.names = FALSE, quote = FALSE)
  bedPeaks <- data.frame(peaks$seqnames, peaks$start, peaks$end)
  write.table(bedPeaks, file = paste0(outDir, 'peaksGroupedBy_', name, '_formatted.bed'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = '\t')
  runBedSort(paste0(outDir, 'peaksGroupedBy_', name, '_formatted.bed'), bedSort)
  runBedToBigBed(paste0(outDir, 'peaksGroupedBy_', name, '_sorted.bed'), bedToBigBed, chromSizes)
}

runBedSort <- function(bedFile, bedSort)
{
  system(paste0(bedSort, ' ', bedFile, ' ', str_replace(bedFile, '_formatted.bed', '_sorted.bed')))
}

runBedToBigBed <- function(bedFile, bedToBigBed, chromSizes)
{
  system(paste0(bedToBigBed, ' ', bedFile, ' ', chromSizes, ' ', str_replace(bedFile, '_sorted.bed', '.bb')))
}

writeGroupPeaksAndBigBed <- function(outDir, peaks, bedSort, bedToBigBed, chromSizes)
{
  lapply(unique(peaks$group), function(x){dir.create(paste0(outDir, x))})
  peaks %>% group_by(group) %>% group_walk(~writePeaksAndRunBigBed(paste0(outDir, .y$group, '/'), ., .y$group, bedSort, bedToBigBed, chromSizes))
}


project <- loadArchRProject(args$Project)

project$celltype <- rep('celltype', times = nCells(project))

project <- addGroupCoverages(project, groupBy = args$groupColName, force = TRUE)

project <- addReproduciblePeakSet(project, groupBy = args$groupColName, pathToMacs2 = findMacs2(), cutOff = args$cutOff, force = TRUE)

project <- addPeakMatrix(project, force = TRUE)

peaks <- getPeaks(project)

dir.create(args$outDir, recursive = TRUE)

writeGroupPeaksAndBigBed(args$outDir, peaks, args$bedSort, args$bedToBigBed, args$chromSizes)

saveArchRProject(project)
