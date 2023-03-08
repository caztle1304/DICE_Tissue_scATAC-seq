library(ArchR)
library(argparse)
library(stringr)

args <- ArgumentParser()
args$add_argument('--Project', default = NULL, type = 'character', help = 'Path to ArchR project')
args$add_argument('--outProjectDir', default = NULL, type = 'character', help = 'Path where new sub projects will be written')
args$add_argument('--outPeakDir', default = NULL, type = 'character', help = 'Path to directory where peaks will be written')
args$add_argument('--groupBy', default = NULL, type = 'character', help = 'Name of column containing each subset to make')
args$add_argument('--bedSort', default = NULL, type = 'character', help = 'Path to bedSort executable')
args$add_argument('--bedToBigBed', default = NULL, type = 'character', help = 'Path to bedToBigBed executable')
args$add_argument('--chromSizes', default = NULL, type = 'character', help = 'Path to chromSizes file for bedToBigBed')
args$add_argument('--minCells', default = 25, type = 'integer', help = 'minCells parameter in addReproduciblePeakSet function for peak calling')

args <- args$parse_args()

subsetProject <- function(project, groupBy, group, outDir) # Subsets an archr project using a reference columna and a group in that column
{ 
  cellColData <- as.data.frame(getCellColData(project, select = groupBy))
  cells2keep <- project$cellNames[which(cellColData[[groupBy]] == group)]
  subsettedProject <- subsetArchRProject(project, cells = cells2keep, outputDirectory = outDir, force = TRUE)
  return(subsettedProject)
}

getPeaks <- function(project) # Gets peaks from a project and formats a data frame to write a bed file
{
  peaks <- getPeakSet(project)
  names(peaks) <- NULL
  peaks <- as.data.frame(peaks)
  return(peaks)
}

runBedSort <- function(bedFile, bedSort) # Runs bedSort from command line
{
  system(paste0(bedSort, ' ', bedFile, ' ', str_replace(bedFile, '_formatted.bed', '_sorted.bed')))
}

runBedToBigBed <- function(bedFile, bedToBigBed, chromSizes) # Runs bedToBigBed from command line 
{
  system(paste0(bedToBigBed, ' ', bedFile, ' ', chromSizes, ' ', str_replace(bedFile, '_sorted.bed', '.bb')))
}

writeGroupPeaksAndBigBed <- function(outDir, peaks, name, bedSort, bedToBigBed, chromSizes) # Writes bed file and runs necessary functions to get final bigBed  
{
  write.csv(peaks, file = paste0(outDir, '/peaksGroupedBy_', name, '_full.csv'), row.names = FALSE, quote = FALSE)
  bedPeaks <- data.frame(peaks$seqnames, peaks$start, peaks$end)
  write.table(bedPeaks, file = paste0(outDir, '/peaksGroupedBy_', name, '_formatted.bed'), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = '\t')
  runBedSort(paste0(outDir, '/peaksGroupedBy_', name, '_formatted.bed'), bedSort)
  runBedToBigBed(paste0(outDir, '/peaksGroupedBy_', name, '_sorted.bed'), bedToBigBed, chromSizes)
}

## main()

project <- loadArchRProject(args$Project)

cellColData <- getCellColData(project)

cellColData <- as.data.frame(cellColData)

for(s in unique(cellColData[[args$groupBy]])) # Iterates each group in the project, creates sub projects and creates bigBed file from each sub project
{
  outProjectDir <- paste0(args$outProjectDir, '/', s)
  dir.create(args$outProjectDir, recursive = TRUE) # Creates necessary sub directories
  outPeakDir <- paste0(args$outPeakDir, '/', s) 
  dir.create(outPeakDir, recursive = TRUE) # Creates necessary sub directories 
  subProject <- subsetProject(project, args$groupBy, s, outProjectDir)
  subProject <- addGroupCoverages(subProject, groupBy = args$groupBy, force = TRUE)
  subProject <- addReproduciblePeakSet(subProject, groupBy = args$groupBy, force = TRUE, minCells = args$minCells)
  subProject <- addPeakMatrix(subProject, force = TRUE)
  saveArchRProject(subProject)
  peaks <- getPeaks(project)
  writeGroupPeaksAndBigBed(outPeakDir, peaks, s, args$bedSort, args$bedToBigBed, args$chromSizes)
}
