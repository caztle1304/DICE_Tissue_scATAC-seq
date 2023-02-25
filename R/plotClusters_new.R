############    - Plotting UMAP   -    ############

# ---
# Authors: Angel Adrian De la Cruz Castillo, Job Rocha
# Date: 2022-10-10
# ---

### -------------------------- Description -------------------------- ###
# Using a UMAP in an ArchR project, this script will plot it colored by different QC metrics: TSS enrichment, number of fragments, sample of origin, reads in promoter, 
# promoter ratio, reads in TSS; as well as by the found clusters.

library(ArchR)
library(argparse)
library(stringr)
set.seed(1)
addArchRThreads(threads = 1)


##Define Arguments####
args<-ArgumentParser()
args$add_argument("--Project",default=NULL,type="character",help="Project path")
args$add_argument("--UMAP_name",default=NULL,type="character",help="UMAP name")
args$add_argument('--Harmony_name', default = NULL, type = 'character', help = 'Harmony correction name')
args$add_argument("--Cluster_name",default=NULL,type="character",help="Clustering name")
args$add_argument("--LSI_name",default=NULL,type="character",help="LSI name")
args$add_argument("--folder",default=NULL,type="character",help="Output path")
args$add_argument('--sampleColName', default = 'Sample', type = 'character', help = 'Name of column where names of samples are stored')
args$add_argument('--resolutions', default = NULL, type = 'character', help = 'Resolutions of clustering to plot separated by commas. Example: 0.1,0.2,0.3,0.4')
args<-args$parse_args()

##Plot Functions###
plot_QCmetrics<-function(df,color_by,colors="viridis"){
  p<-ggplot(data=df,aes(x=.data[["Dimension_1"]],y=.data[["Dimension_2"]],color=.data[[color_by]]))+geom_point()+scale_color_continuous(type=colors)+theme_bw()+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  return(p)
}

resolutions <- str_split(args$resolutions, ',')[[1]]

##Main###
#Retrieve Data regular UMAP#
Project<-loadArchRProject(args$Project)
metrics<-getCellColData(Project)
emb<-getEmbedding(Project,embedding=args$UMAP_name)
cells<-intersect(rownames(metrics),rownames(emb))

metrics<-metrics[order(cells),]
emb<-emb[order(cells),]
metrics["Dimension_1"]<-emb[,paste0(args$LSI_name,"#UMAP_Dimension_1")]
metrics["Dimension_2"]<-emb[,paste0(args$LSI_name,"#UMAP_Dimension_2")]
metrics<-as.data.frame(metrics)

##Plots##
samples <- plotEmbedding(ArchRProj = Project, colorBy = "cellColData", name = args$sampleColName, embedding = args$UMAP_name)

mfrags<-metrics[metrics[,"nFrags"]<=quantile(metrics[,"nFrags"],0.95),]
pfrags<-plot_QCmetrics(mfrags,"nFrags","viridis")

mTSS<-metrics[metrics[,"TSSEnrichment"]<=quantile(metrics[,"TSSEnrichment"],0.95),]
pTSS<-plot_QCmetrics(mTSS,"TSSEnrichment","viridis")

mreads.TSS<-metrics[metrics[,"ReadsInTSS"]<=quantile(metrics[,"ReadsInTSS"],0.95),]
preads.TSS<-plot_QCmetrics(mreads.TSS,"ReadsInTSS","viridis")

mreads.Prom<-metrics[metrics[,"ReadsInPromoter"]<=quantile(metrics[,"ReadsInPromoter"],0.95),]
preads.Prom<-plot_QCmetrics(mreads.Prom,"ReadsInPromoter","viridis")

mPR<-metrics[metrics[,"PromoterRatio"]<=quantile(metrics[,"PromoterRatio"],0.95),]
pPR<-plot_QCmetrics(mPR,"PromoterRatio","viridis")

mNR<-metrics[metrics[,"NucleosomeRatio"]<=quantile(metrics[,"NucleosomeRatio"],0.95),]
pNR<-plot_QCmetrics(mNR,"NucleosomeRatio","viridis")

## Regular UMAP plot

dir.create(args$folder,recursive=TRUE)
pdf(paste0(args$folder,"/Clustering_", args$Cluster_name, '.pdf'))
samples
for (res in resolutions)
{
  print(plotEmbedding(ArchRProj = Project, colorBy = 'cellColData', name = paste0(args$Cluster_name, '_', res), embedding = args$UMAP_name, reducedDims = args$LSI_name ))

}
pfrags
pTSS
preads.TSS
preads.Prom
pPR
pNR
dev.off()



#Retrieve Data Harmony UMAP#
emb<-getEmbedding(Project,embedding=paste0(args$Harmony_name, '_UMAP'))
cells<-intersect(rownames(metrics),rownames(emb))

metrics<-metrics[order(cells),]
emb<-emb[order(cells),]
metrics["Dimension_1"]<-emb[,paste0(args$Harmony_name,"#UMAP_Dimension_1")]
metrics["Dimension_2"]<-emb[,paste0(args$Harmony_name,"#UMAP_Dimension_2")]
metrics<-as.data.frame(metrics)

##Plots##
samples <- plotEmbedding(ArchRProj = Project, colorBy = "cellColData", name = args$sampleColName, embedding = paste0(args$Harmony_name, '_UMAP'))

mfrags<-metrics[metrics[,"nFrags"]<=quantile(metrics[,"nFrags"],0.95),]
pfrags<-plot_QCmetrics(mfrags,"nFrags","viridis")

mTSS<-metrics[metrics[,"TSSEnrichment"]<=quantile(metrics[,"TSSEnrichment"],0.95),]
pTSS<-plot_QCmetrics(mTSS,"TSSEnrichment","viridis")

mreads.TSS<-metrics[metrics[,"ReadsInTSS"]<=quantile(metrics[,"ReadsInTSS"],0.95),]
preads.TSS<-plot_QCmetrics(mreads.TSS,"ReadsInTSS","viridis")

mreads.Prom<-metrics[metrics[,"ReadsInPromoter"]<=quantile(metrics[,"ReadsInPromoter"],0.95),]
preads.Prom<-plot_QCmetrics(mreads.Prom,"ReadsInPromoter","viridis")

mPR<-metrics[metrics[,"PromoterRatio"]<=quantile(metrics[,"PromoterRatio"],0.95),]
pPR<-plot_QCmetrics(mPR,"PromoterRatio","viridis")

mNR<-metrics[metrics[,"NucleosomeRatio"]<=quantile(metrics[,"NucleosomeRatio"],0.95),]
pNR<-plot_QCmetrics(mNR,"NucleosomeRatio","viridis")

### Harmony-corrected UMAP

pdf(paste0(args$folder,"/Clustering_",args$Harmony_name, '.pdf'))
samples
for (res in resolutions)
{
  print(plotEmbedding(ArchRProj = Project, colorBy = 'cellColData', name =paste0(args$Harmony_name, '_', res), embedding = paste0(args$Harmony_name, '_UMAP'), reducedDims = args$Harmony_name))

}
pfrags
pTSS
preads.TSS
preads.Prom
pPR
pNR

dev.off()
