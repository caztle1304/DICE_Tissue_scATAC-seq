############    -  Processing of scATAC-seq data  -    ############

# ---
# Author: Angel Adrian De la Cruz Castillo
# Date: 2022-10-14
# ---

### -------------------------- Description -------------------------- ###
# Based on the result(s) of the clustering(s), this program will plot: cells per cluster, proportion of each cluster, proportion of each sample in ArchR project in a given cluster,
# total proportion of each sample in project, and will create a table with ArchR QCs grouped by sample

library(ArchR)
library(argparse)
library(stringr)
library(dplyr)
library(ggrepel)


args <- ArgumentParser()
args$add_argument('--Project', default = NULL, type = 'character', help = 'Project path')
args$add_argument('--plotDir', default = NULL, type = 'character', help = 'Path to directory to contain output plots')
args$add_argument('--clusterName', default = NULL, type = 'character', help = 'Name of clustering')
args$add_argument('--Harmony_name', default = 'Harmony', type = 'character', help = 'Name of harmony dim correction')
args$add_argument('--resolutions', default = NULL, type = 'character', help = 'Resolutions in clustering to use')
args$add_argument('--sampleColName', default = 'Sample', type = 'character', help = 'Column of project where sample data is stored')
args <- args$parse_args()



plotCellsPerCluster <- function(project, clusterName, celltype)
{
  cellColData <- getCellColData(project)
  cellColData <- as.data.frame(cellColData)
  clusters <- as.data.frame(table(cellColData[, get('clusterName')]))
  clusters <- clusters[str_order(clusters$Var1, numeric = TRUE),]
  clusters$Var1 <- factor(clusters$Var1, levels = clusters$Var1)

  plot <- ggplot(data = clusters, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = Freq), vjust = 1.6, color = 'black', size = 3.5) +
  ggtitle(paste0(clusterName, ', nCells = ', sum(clusters$Freq))) +
  theme(legend.position = 'none') + xlab('Cluster') + ylab('nCells')

  return(plot)
}

plotSampleProportionsPerCluster <- function(project, clusterName, sampleColName)
{
  data <- getCellColData(project, select = c(get('sampleColName'), get('clusterName')))
  data <- as.data.frame(data)
  df <- data.frame(data %>% group_by(get(get('sampleColName'))) %>% count(get(get('clusterName'))))
  setnames(df, new = c('Sample', 'Cluster', 'N'))

  plot <- ggplot(df, aes(fill = Sample, y = Cluster, x = N)) + geom_bar(position = 'fill', stat = 'identity') +
  coord_flip() + ggtitle(paste0(clusterName, ' Clustering Sample Proportions'))

  return(plot)
}


getArchRSummary <- function(project, sampleCol){
  data <- getCellColData(project)
  data <- as.data.frame(data)
  summary <- data %>% group_by(get(get('sampleCol'))) %>% summarise_at(c('nMultiFrags', 'nMonoFrags', 'nFrags', 'nDiFrags', 'TSSEnrichment', 'ReadsInTSS', 'ReadsInPromoter', 'ReadsInBlacklist', 'PromoterRatio', 'NucleosomeRatio', 'BlacklistRatio'), mean)
  summary$nCells <- as.data.frame(table(data[,get('sampleCol')]))$Freq
  summary <- cbind(summary[,1], summary$nCells, summary[,-c(1,13)])
  setnames(summary, new = c('Sample', 'nCells', 'nMultiFrags_mean', 'nMonoFrags_mean', 'nFrags_mean', 'nDiFrags_mean', 'TSSEnrichment_mean', 'ReadsInTSS_mean', 'ReadsInPromoter_mean', 'ReadsInBlacklist_mean', 'PromoterRatio_mean', 'NucleosomeRatio_mean', 'BlacklistRatio_mean'))

  return(as.data.frame(summary))
}

clusterProportions <- function(project, clusterName)
{
  cellColData <- getCellColData(project)
  data <- as.data.frame(table(cellColData[, get('clusterName')]))
  data <- data[order(data$Freq),]
  data$Var1 <- factor(data$Var1, levels = data$Var1)
  data <- data %>% mutate(Perc = paste0(round(Freq/sum(Freq)*100, digits = 2), '%'), text_y = sum(Freq) - (cumsum(Freq) - Freq/2))
  #csum = rev(cumsum(rev(Freq))),
  #       pos = Freq/2 + lead(csum, 1),
  #       pos = if_else(is.na(pos), Freq/2, pos))
  plot <- ggplot(data, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + # remove background, grid, numeric labels
  ggtitle(clusterName, ' Cluster proportions') + 
  geom_label_repel(aes(label = Perc, y = text_y), nudge_x = .6, nudge_y = .6, show.legend = F, size = 5) +
  guides(fill = guide_legend(title = 'Cluster'))
  return(plot)
}


sampleProportions <- function(project, sampleColName)
{
  cellColData <- getCellColData(project)
  data <- as.data.frame(table(cellColData[, get('sampleColName')]))
  data <- data [order(data$Freq),]
  data$Var1 <- factor(data$Var1, levels = data$Var1)
  data <- data %>% mutate(Perc = paste0(round(Freq/sum(Freq)*100, digits = 2), '%'), text_y = sum(Freq) - (cumsum(Freq) - Freq/2))
  plot <- ggplot(data, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  ggtitle('Total sample proportions') +
  geom_label_repel(aes(label = Perc, y = text_y), nudge_x = .6, nudge_y = .6, show.legend = F, size = 5) +
  guides(fill = guide_legend(title = 'Sample'))
  return(plot)
}

resolutions <- str_split(args$resolutions, ',')[[1]]
project <- loadArchRProject(args$Project)


## Regular UMAP clustering

pdf(paste0(args$plotDir, 'clusterStats.pdf'))

for (res in resolutions)
{
  clustering <- paste0(args$clusterName, '_', res)
  print(plotCellsPerCluster(project, clustering, args$celltype))
  print(plotSampleProportionsPerCluster(project, clustering, args$sampleColName))
  print(clusterProportions(project, clustering))
}

print(sampleProportions(project, args$sampleColName))

dev.off()


## Harmony UMAP

pdf(paste0(args$plotDir, 'clusterStats_Harmony.pdf'))

# If multiple resolutions provided, it will create plots for each one of them
for (res in resolutions)
{
  clustering <- paste0(args$Harmony_name, '_', res)
  print(plotCellsPerCluster(project, clustering, args$celltype))
  print(plotSampleProportionsPerCluster(project, clustering, args$sampleColName))
  print(clusterProportions(project, clustering))
}

print(sampleProportions(project, args$sampleColName))

dev.off()


summary <- getArchRSummary(project, args$sampleColName)

write.csv(summary, file = paste0(args$plotDir, 'projectSummary.csv'), quote = FALSE, row.names = FALSE)

