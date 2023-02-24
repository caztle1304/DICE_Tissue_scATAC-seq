library(ArchR)
library(stringr)
library(rjson)
library(argparse)
set.seed(1)

args<-ArgumentParser()
args$add_argument("--Project",default=NULL,type="character",help="Project path")
args$add_argument("--useMatrix",default="TileMatrix",type="character",help="Matrix to use for clustering. Default: TileMatrix")
args$add_argument("--IterativeLSI_name",default=NULL,type="character",help="Name for LSI")
args$add_argument("--iterations",default=2,type="integer",help="Number of iterations. Default: 2")
args$add_argument("--LSI_resolution",default=0.2,type="double",help="Resolution for LSI. Default: 0.2")
args$add_argument("--sampleCells",default=10000,type="integer",help="Number of cells to sample. Default: 10,000")
args$add_argument("--nstart",default=10,type="integer",help="Default: 10")
args$add_argument("--varFeatures",default=25000,type="integer",help="Number of variable features to use. Default: 25000")
args$add_argument("--dimsToUse",default="1_30",type="character",help="Dimensions to use. Default: 1_30")
args$add_argument("--method",default="Seurat",type="character",help="Method to use for clustering. Default: Seurat")
args$add_argument("--clustering_name",default=NULL,type="character",help="Name for current clustering analysis")
args$add_argument("--cluster_resolutions",default='0.1,0.2,0.4,0.6,0.8',type="character",help="Cluster Resolution. Default: 0.1,0.2,0.4,0.6,0.8")
args$add_argument("--UMAP_name",default=NULL,type="character",help="Name for UMAP")
args$add_argument("--nNeighbors",default=30,type="integer",help="Number of neighbours in UMAP. Default: 30")
args$add_argument("--minDist",default=0.5,type="double",help="min distance in UMAP. Default: 0.5")
args$add_argument("--metric",default="cosine",type="character",help="UMAP metric. Default: cosine")
args$add_argument("--threads",default=5,type="integer",help="number of threads. Default: 5")
args$add_argument('--Harmony_name', default = 'Harmony', type = 'character', help = 'Name of harmony dim reduction')
args$add_argument('--groupby', default = NULL, type = 'character', help = 'Group by parameter in harmony, i.e. samples in project')

args<-args$parse_args()

addArchRThreads(threads = 1)

Project<-loadArchRProject(args$Project)
dims<-args$dimsToUse
dims<-str_split(dims,"_")
dims<-as.integer(dims[[1]][1]):as.integer(dims[[1]][2])

resolutions <- str_split(args$cluster_resolutions, ',')[[1]]

Project <- addIterativeLSI(
    ArchRProj = Project,
    useMatrix = args$useMatrix,
    name = args$IterativeLSI_name,
    iterations = args$iterations,
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(args$LSI_resolution),
        sampleCells = args$sampleCells,
        n.start = args$nstart
    ),
    varFeatures = args$varFeatures,
    dimsToUse = dims,
    force=TRUE
)
Project <- addHarmony(ArchRProj = Project,
reducedDims = args$IterativeLSI_name,
name = args$Harmony_name,
groupBy = args$groupby, 
force = TRUE)


for (res in resolutions){

Project <- addClusters(
    input = Project,
    reducedDims = args$IterativeLSI_name,
    method = args$method,
    name = paste0(args$clustering_name, '_', res),
    resolution = as.double(res),
    force=TRUE)

Project <- addClusters(input = Project,
reducedDims = args$Harmony_name,
method = args$method,
name = paste0(args$Harmony_name, '_', res),
resolution = as.double(res),
force = TRUE)
}


Project <- addUMAP(
    ArchRProj = Project,
    reducedDims = args$IterativeLSI_name,
    name = args$UMAP_name,
    nNeighbors = args$nNeighbors,
    minDist = args$minDist,
    metric = args$metric,
    force=TRUE
)

Project <- addUMAP(
    ArchRProj = Project,
    reducedDims = args$Harmony_name,
    name = paste0(args$Harmony_name, '_UMAP'),
    nNeighbors = args$nNeighbors,
    minDist = args$minDist,
    metric = args$metric,
    force=TRUE
)


data<-data.frame(arguments=names(args),values=unlist(args))

dir.create(paste0(args$Project,"/ClusterParams/"),recursive=TRUE)
paramfile<-paste0(args$Project,"/ClusterParams/",args$clustering_name,".csv")
write.table(data,paramfile,row.names=FALSE,sep=',')

saveArchRProject(ArchRProj = Project, outputDirectory = args$Project, load = FALSE)
