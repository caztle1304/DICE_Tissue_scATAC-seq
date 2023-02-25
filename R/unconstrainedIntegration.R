############    - scATAC-seq integration with scRNA-seq  -    ############

# ---
# Author:Angel Adrian De la Cruz Castillo
# Date: 2022-12-16
# ---

### -------------------------- Description -------------------------- ###
# Using an ArchR project and a Seurat object containing scRNA-seq data, this script will perform unconstrained integration and will plot the results 
# in a UMAP, as well as the prediction score and some plots for QC: violoin plot with prediction score per cell type, number and proportion of each ell type

library(plyr)
library(ArchR)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(Seurat)
library(pheatmap)
library(argparse)
library(stringr)
library(dplyr)
library(ggrepel)
set.seed(1)

args <- ArgumentParser()

args$add_argument('--Project', default = NULL, type = 'character', help = 'Path to archr project ')
args$add_argument('--SeuratObject', default = NULL, type = 'character', help = 'Path to Seurat object to integrate')
args$add_argument('--SoClustering', default = NULL, type = 'character', help = 'Name of column in seurat object containing clusters')
args$add_argument('--SoIdentities', default = NULL, type = 'character', help = 'Name of identities corresponding to each cluster in order and separated by commas, that means, the first ID has to match with cluster 0, the second ID with cluster 1 and so on')
args$add_argument('--ReducedDims', default = NULL, type = 'character', help = 'Name of LSI/Harmony in ArchR project')
args$add_argument('--clusterName', default = NULL, type = 'character', help = 'Base name of clustering in ArchR project')
args$add_argument('--UMAP', default = NULL, type = 'character', help = 'UMAP in ArchR project')
args$add_argument('--resolutions', default = NULL, type = 'character', help = 'Comma-separated values of resolutions of atac clustering for confusion matrices')
args$add_argument('--plotDir', default = NULL, type = 'character', help = 'Plot directory')
args$add_argument('--PredictionSuffix', default = NULL, type = 'character', help = 'Suffix to add to new prediction score, group and cell columns in ArchR project')
args$add_argument('--ptSize', default = 0.1, type = 'double', help = 'Size of points in UMAP')
args$add_argument('--nCells_ATAC', default = 10000, type = 'integer', help = 'Number of cells to sample from ATAC to perform integration')
args$add_argument('--nCells_RNA', default = 10000, type = 'integer', help = 'number of cells to sample from RNA to perform integration')
args$add_argument('--plotAs', default = 'points', type = 'character', help = '"hex" or "points", if points dot plot, if hex, hexaplot')

args <- args$parse_args()


changeSeuratGeneNames <- function(so)
{ ## Changes all possible Ensembl IDs to gene symbols
  genesDict <- ensembldb::select(EnsDb.Hsapiens.v86, keys = rownames(so), keytype = 'GENEID', columns = 'SYMBOL')
  newGenes <- mapvalues(rownames(so), from = genesDict$GENEID, to = genesDict$SYMBOL)
  # Re-names necessary slots
  rownames(so@assays$RNA@counts) <- newGenes
  rownames(so@assays$RNA@data) <- newGenes
  return(so)
}


removeSeuratGenes <- function(so, archrproj)
{ ## Removes genes from seurat object to keep only those also present in archr
  counts <- GetAssayData(so, assay = 'RNA')
  atac_features <- getFeatures(archrproj)
  counts <- counts[rownames(so) %in% atac_features,]
  so <- subset(so, features = rownames(counts))
  ## Adds rownames to meta.features to avoid error
  so[['RNA']]@meta.features <- data.frame(row.names = rownames(so[['RNA']]))
  return(so)
}



## main()

so_rna <- readRDS(args$SeuratObject)
atac <- loadArchRProject(args$Project)

## Adding additional column to seurat object with cell identities
so_rna$cellIDs <- mapvalues(as.vector(so_rna[[args$SoClustering]][[args$SoClustering]]), from = sort(unique(so_rna[[args$SoClustering]])[[args$SoClustering]]), to = unlist(str_split(args$SoIdentities, ',')))

so_rna <- changeSeuratGeneNames(so_rna)
so_rna <- removeSeuratGenes(so_rna, atac)

atac <- addGeneIntegrationMatrix(atac, useMatrix = 'GeneScoreMatrix',
reducedDims = args$ReducedDims,
seRNA = so_rna,
addToArrow = TRUE,
groupRNA = 'cellIDs',
force = TRUE,
sampleCellsATAC = args$nCells_ATAC,
sampleCellsRNA = args$nCells_RNA,
nameCell = paste0('predictedCell_', args$PredictionSuffix),
nameGroup = paste0('predictedGroup_', args$PredictionSuffix),
nameScore = paste0('predictedScore_', args$PredictionSuffix))

ifelse(args$plotAs == 'points', ptSize <- args$ptSize, ptSize <- NULL)

emb1 <- plotEmbedding(atac, embedding = args$UMAP, colorBy = 'cellColData', name = paste0('predictedGroup_', args$PredictionSuffix), size = ptSize, plotAs = args$plotAs)
emb2 <- plotEmbedding(atac, embedding = args$UMAP, colorBy = 'cellColData', name = paste0('predictedScore_', args$PredictionSuffix), size = ptSize, plotAs = args$plotAs)

dir.create(args$plotDir, recursive = TRUE)

pdf(paste0(args$plotDir, 'UnconstrainedIntegration.pdf'))

print(emb1)
print(emb2)

dev.off()

saveArchRProject(atac)
