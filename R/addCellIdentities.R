############    -  Add annotations of subpopulations to ArchR project  -    ############

# ---
# Author: Angel Adrian De la Cruz Castillo
# Date: 2022-11-14
# ---

### -------------------------- Description -------------------------- ###
# Based on the chosen resolution of clustering, this program will add an additional column of data to an ArchR project specifying the subpopulation to which 
# a cell belongs

library(ArchR)
library(argparse)
library(stringr)
library(plyr)

args <- ArgumentParser()
args$add_argument('--clusters', default = NULL, type = 'character', help = 'Clusters to add identity to; must be comma-separated and match the order of cell identities')
args$add_argument('--celltypes', default = NULL, type = 'character', help = 'Cell identities to add to clusters, must match the order of clusters and be comma-separated')
args$add_argument('--Project', default = NULL, type = 'character', help = 'Path to project')
args$add_argument('--clusterName', default = NULL, type = 'character', help = 'Name of clustering in project to use')
args$add_argument('--newColName', default = NULL, type = 'character', help = 'Name of new column to be created when adding cell identities ')
args <- args$parse_args()

project <- loadArchRProject(args$Project)

clusters <- str_split(args$clusters, ',')[[1]]

celltypes <- str_split(args$celltypes, ',')[[1]]

cellColData <- getCellColData(project)

cellColData <-as.data.frame(cellColData)

cellIdentities <- mapvalues(cellColData[, args$clusterName], from = clusters, to =celltypes)

project <- addCellColData(project, data = cellIdentities, name = args$newColName, cells = project$cellNames, force = TRUE)

saveArchRProject(project)
