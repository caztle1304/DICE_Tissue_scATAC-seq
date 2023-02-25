############    -  Gene score calculation and visualization in UMAP  -    ############

# ---
# Author: Angel Adrian De la Cruz Castillo
# Date: 2022-12-15
# ---

### -------------------------- Description -------------------------- ###
# This script will calculate gene scores of an ArchR project, impute them, and will visualize them on regular and Harmony UMAPs using a list of user-provided marker genes.

library(ArchR)
library(argparse)
library(stringr)
set.seed(1)
addArchRThreads(threads = 1)

args<-ArgumentParser()
args$add_argument("--Project",default=NULL,type="character",help="Project path")
args$add_argument("--UMAP_name",default=NULL,type="character",help="UMAP name")
args$add_argument("--LSI_name",default=NULL,type="character",help="LSI name")
args$add_argument("--folder",default=NULL,type="character",help="Output path")
args$add_argument("--plot_name",default=NULL,type="character",help="Plot name")
args$add_argument("--genes",default="",type="character",help="Marker genes for GeneScore plots")
args$add_argument('--Harmony_name', default = NULL, type = 'character', help = 'name of harmony dim reduction')
args$add_argument('--plotAs', default = 'hex', type = 'character', help = '"points" for dots or "hex" for hexaplot, default:"hex"')
args$add_argument('--size', default = 0.1, type = 'double', help = 'If plotAs is set to "points", defines size of the dots, Default: 0.1')
args<-args$parse_args()


plotAs <- args$plotAs
ifelse(plotAs == 'points', size <- args$size, size <- NULL)

Project<-loadArchRProject(args$Project)
marker_genes<-scan(args$genes,what="list",sep="\n")

## Regular UMAP

MGnoimp <- plotEmbedding(
  ArchRProj = Project,
  colorBy = "GeneScoreMatrix",
  name = marker_genes,
  embedding = args$UMAP_name, 
  plotAs = plotAs, 
  size = size)

Project <- addImputeWeights(Project,reducedDims = args$LSI_name)

MGimp<-plotEmbedding(
  ArchRProj = Project,
  colorBy = "GeneScoreMatrix",
  name = marker_genes,
  embedding = args$UMAP_name,
  imputeWeights = getImputeWeights(Project), 
  plotAs = plotAs, 
  size = size
)


dir.create(args$folder,recursive=TRUE)

pdf(paste0(args$folder,"/",args$plot_name,"_withimputation.pdf"))
MGimp
dev.off()

pdf(paste0(args$folder,"/",args$plot_name,"_noimputation.pdf"))
MGnoimp
dev.off()

## Harmony UMAP 

 MGnoimp_harmony <- plotEmbedding(ArchRProj = Project, 
   colorBy = 'GeneScoreMatrix', 
   name = marker_genes,
 embedding = paste0(args$Harmony_name, '_UMAP'), 
   plotAs = plotAs, 
   size = size)

Project <- addImputeWeights(Project, reducedDims = args$Harmony_name)

 MGimp_harmony <- plotEmbedding(ArchRProj = Project,
 colorBy = 'GeneScoreMatrix', 
name = marker_genes,
embedding = paste0(args$Harmony_name, '_UMAP'),
imputeWeights = getImputeWeights(Project),
plotAs = plotAs, 
size = size)


pdf(paste0(args$folder, '/', args$plot_name, '_withimputation_Harmony.pdf'))
MGimp_harmony
dev.off()

pdf(paste0(args$folder, '/', args$plot_name, '_noimputation_Harmony.pdf'))
MGnoimp_harmony
dev.off()


saveArchRProject(Project)
