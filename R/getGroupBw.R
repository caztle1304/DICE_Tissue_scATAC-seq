############    - BigWig file creation from ArchR project  -    ############

# ---
# Author: Angel Adrian De la Cruz Castillo
# Date: 2023-02-10
# ---

### -------------------------- Description -------------------------- ###
# Based on a given column in an ArchR project, this script will create bigWig files for each group, move them to an output directory and create 
# necessary subdirectories in the process

library(ArchR)
library(dplyr)
library(argparse)
library(stringr)

args <- ArgumentParser()
args$add_argument('--Project', default = NULL, type = 'character', help = 'Path to ArchR project')
args$add_argument('--outDir', default = NULL, type = 'character', help = 'Path to directory to save bw files')
args$add_argument('--groupBy', default = NULL, type = 'character', help = 'Column in project that specifies how to group bigwig files. If celltype, will create a new temporary column to generate pseudo-bulk based on all cells in project. If multiple columns, must be comma-separated')
args$add_argument('--normMethod', default = NULL, type = 'character', help = 'Normalization method to use for bw file, if multiple must be comma-separated, Accepted: ReadsInTSS, ReadsInPromoter, nFrags, nCells')
args$add_argument('--tileSize', default = NULL, type = 'character', help = 'Size of tile to account for accesibility in a given genomic region; if multiple, must be comma-separated')


args <- args$parse_args()

getBigWig <- function(project, normMethod, outDir, groupBy, tileSize)
{
  getGroupBW(project, normMethod = normMethod, tileSize = tileSize, groupBy = groupBy, maxCells = NULL)
  inputFile <- paste0(getOutputDirectory(project), '/GroupBigWigs/', groupBy, '/*.bw')
  dir.create(outDir, recursive = TRUE)
  system(paste0('mv ', inputFile, ' ', outDir))
}

## main()

project <- loadArchRProject(args$Project)


project$celltype <- rep('celltype', times = nCells(project))


groupBy <- unlist(str_split(args$groupBy, ','))
norm <- unlist(str_split(args$normMethod, ','))
tileSize <- unlist(str_split(args$tileSize, ','))

mainOutDir <- args$outDir

for (g in groupBy)
{
  for(n in norm)
  {
    for(ts in tileSize)
    { 
      outDir <- paste0(mainOutDir, '/by_', g, '/normBy_', n, '/TileSize_', ts, '/')
      getBigWig(project, n, outDir, g, as.numeric(ts))
    }
  }
}
