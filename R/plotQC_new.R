############    -  Plotting QC metrics  -    ############

# ---
# Author: Angel Adrian De la Cruz Castillo
# Date: 2022-10-10
# ---

### -------------------------- Description -------------------------- ###
#This script will plot TSS enrichment and distribution of fragment size, as well as a TSS vs Fragments plot of an ArchR project
library(ArchR)
library(argparse)
library(stringr)
set.seed(1)
addArchRThreads(threads = 10)

args<-ArgumentParser()
args$add_argument("--Project",default="",type="character",help="Absolute path to Project")
args$add_argument("--samples",default="all",type="character",help="Samples to use")
args$add_argument("--TSS",default=4,type="integer",help="min TSS value to filter cells")
args$add_argument("--Frags",default=1000,type="integer",help="min number of Frags to filter cells")
args$add_argument("--groupby", default="Sample",type="character",help="Group by parameter in plot")
args$add_argument("--folder",default=NULL,type="character",help="Output directory name")
args$add_argument("--sampleColName", default = 'Sample', type = 'character', help = 'Name of slot containing sample info')
args$add_argument("--name", default  = NULL, type = 'character', help = 'Name of the prefix of plots, i.e. name of your project')
args<-args$parse_args()


samples<-str_split(args$samples,",")[[1]]



Project<-loadArchRProject(args$Project)
df <- getCellColData(Project, select = c("log10(nFrags)", "TSSEnrichment",args$sampleColName))
if (!args$samples=="all"){
        df <- df[df$Sample %in% samples,c("log10(nFrags)", "TSSEnrichment")]}


TSSvsFrags <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = args$TSS, lty = "dashed") + geom_vline(xintercept = log10(args$Frags), lty = "dashed")



pdf(paste0(args$folder,"/",args$name,"_TSSvsFrags.pdf"))
TSSvsFrags
dev.off()

pdf(paste0(args$folder,"/",args$name,"_Fragment_Size_Distribution.pdf"))
plotFragmentSizes(ArchRProj = Project, groupBy = args$groupby)
dev.off()

pdf(paste0(args$folder,"/",args$name,"_TSSEnrichment.pdf"))
plotTSSEnrichment(ArchRProj = Project, args$groupby, threads = 10)
dev.off()
