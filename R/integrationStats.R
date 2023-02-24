library(ArchR)
library(ggrepel)
library(pheatmap)
library(stringr)
library(argparse)
library(dplyr)

args <- ArgumentParser()
args$add_argument('--Project', default = NULL, type = 'character', help = 'Path to ArchR project')
args$add_argument('--resolutions', default = NULL, type = 'character', help = 'Resolutions for confusion matrix')
args$add_argument('--plotDir', default = NULL, type = 'character', help = 'Path to directory of output plots')
args$add_argument('--suffix', default = NULL, type = 'character', help = 'Suffix of prediciton columns')
args$add_argument('--clusterName', default = NULL, type = 'character', help = 'Base name of clustering to use in confusion matrix')
args <- args$parse_args()


plotPredictionScorePerCelltype <- function(archrproj, suffix)
{
  ## Plots violin plot with boxplot for scoring of each cell type
cellColData <- getCellColData(archrproj, select = paste0(c('predictedGroup_', 'predictedScore_'), suffix))
cellColData <- as.data.frame(cellColData)

p <- ggplot(cellColData, aes(x = get(paste0('predictedGroup_', suffix)), y = get(paste0('predictedScore_', suffix)), fill = get(paste0('predictedGroup_', suffix)))) +
     geom_violin(width = 0.5, alpha = 0.2) +
     geom_boxplot(width = 0.1) +
     theme_bw() +
     ggtitle('Prediction scores by group') +
     theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
     xlab('Cell type') + ylab('Prediction score')
return(p)
}

plotCellsPerCellType <- function(archrproj, suffix)
{ ## Plots number of cells of each identified cell type
  cellColData <- getCellColData(archrproj, select = paste0('predictedGroup_', suffix))
  cellColData <- as.data.frame(cellColData)
  df <- as.data.frame(table(cellColData))
  plot <- ggplot(df, aes(x = cellColData, y = Freq, fill = cellColData)) + geom_bar(stat = 'identity') +
          geom_text(aes(label = Freq), vjust = 1.6, color = 'black', size = 3.5) + ggtitle('Number of cells per cell type') +
          theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + xlab('Cell Type') + ylab('nCells')
  return(plot)

}

plotCellTypeProportions <- function(archrproj, suffix)
{ ## Plots pie chart showing proportions of each identified cell type
  cellColData <- getCellColData(archrproj, select = paste0('predictedGroup_', suffix))
  df <- as.data.frame(table(as.data.frame(cellColData)))
  df <- df[order(df$Freq),]
  df$Var1 <- factor(df$Var1, levels = df$Var1)
  df <- df %>% mutate(Perc = paste0(round(Freq/sum(Freq)*100, digits = 2), '%'), text_y = sum(Freq) - (cumsum(Freq) - Freq/2)) ## text_y is for position of each label in pie chart
  plot <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + # remove background, grid, numeric labels
  ggtitle('Cell Type proportions') +
  geom_label_repel(aes(label = Perc, y = text_y), nudge_x = .6, nudge_y = .6, show.legend = F, size = 5) +
  guides(fill = guide_legend(title = 'Cell type'))
  return(plot)
}


plotConfusionMatrix <- function(archrproj, clusterName, predictedGroupColumn, title)
{
  cellColData <- getCellColData(archrproj)
  cM <- as.matrix(confusionMatrix(cellColData[[clusterName]], cellColData[[predictedGroupColumn]]))
  cM <- cM / Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(mat = as.matrix(cM),
                          color = paletteContinuous('whiteBlue'),
                          border_color = 'black', main = title)
  return(p)
}


# main()


project <- loadArchRProject(args$Project)

resolutions <- unlist(str_split(args$resolutions, ','))

dir.create(args$plotDir)

pdf(paste0(args$plotDir, 'integrationStats.pdf'))

for(res in resolutions)
{
  print(plotConfusionMatrix(project, paste0(args$clusterName, '_' ,res), paste0('predictedGroup_', args$suffix), paste0('Confusion matrix \n Resolution = ', res)))
}

print(plotPredictionScorePerCelltype(project, args$suffix))

print(plotCellsPerCellType(project, args$suffix))

print(plotCellTypeProportions(project, args$suffix))

dev.off()
