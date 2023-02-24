Rscript=/mnt/BioApps/R/4.0.1/bin/Rscript

#### LSI, Harmony, Clustering, UMAP
Clustering=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/Clustering_new.R
project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k
varFeatures=30000
UMAP=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_UMAP
dims2use=1_40
LSI=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_LSI
resolutions=0.1,0.2,0.4,0.6,0.8
clusterName=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_clustering
harmony=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony
sampleColName=Sample

$Rscript $Clustering --Project $project --varFeatures $varFeatures --UMAP_name $UMAP --dimsToUse $dims2use --IterativeLSI_name $LSI --cluster_resolutions $resolutions --clustering_name $clusterName --Harmony_name $harmony --groupby $sampleColName

#### Plot UMAP colored by clustering and QCs
PlotClusters=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/plotClusters_new.R
sampleColName=Sample
project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k
plotDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/myeloid/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k/
UMAP=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_UMAP
clusterName=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_clustering
LSI=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_LSI
resolutions=0.1,0.2,0.4,0.6,0.8
harmony=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony

$Rscript $PlotClusters --sampleColName $sampleColName --folder $plotDir --Project $project --UMAP_name $UMAP --Cluster_name $clusterName --LSI_name $LSI --resolutions $resolutions --Harmony_name $harmony

#### Plot UMAP colored by gene scores
PlotGeneScores=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/plotGeneScores_new.R
project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k
UMAP=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_UMAP
LSI=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_LSI
plotDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/myeloid/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k/
markerGenes=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/markerGenes/myeloidMarkerGenes.txt
harmony=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony

$Rscript $PlotGeneScores --Project $project --UMAP_name $UMAP --LSI_name $LSI --folder $plotDir --genes $markerGenes --plot_name myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_geneScores.pdf --Harmony_name $harmony

#### Plot general QCs
PlotQC=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/plotQC_new.R
project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k
sampleColName=Sample
plotDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/myeloid/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k/
name=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k

$Rscript $PlotQC --Project $project --TSS 0 --Frags 0 --groupby $sampleColName --folder $plotDir --name $name

#### General and useful statistics
ArchRStats=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/ArchRStats.R
project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k
plotDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/myeloid/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k/
clusterName=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_clustering
harmony=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony
resolutions=0.1,0.2,0.4,0.6,0.8
sampleColName=Sample

$Rscript $ArchRStats --Project $project --plotDir $plotDir --clusterName $clusterName --Harmony_name $harmony --resolutions $resolutions --sampleColName $sampleColName
#### Integration with scRNA-seq
integration=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/unconstrainedIntegration.R
project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k
seuratObject=/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/paper_items/SeuratObj_Healthy-Myeloid.RDS
soIdentities=Classical_monocytes,Non_Classical_monocytes,Tissue_resident_Macrophages,Macrophage_cDC2_mixture,cDC2s,pDCs,cDC1s,mregDC
reducedDims=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony
umap=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_UMAP
resolutions=0.1,0.2,0.4,0.6,0.8
plotDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/integration_roundTwo/myeloid/
PredictionSuffix=Un_Round2
ptSize=0.1
nCells_ATAC=7701
nCells_RNA=100000 ## Originally there were 132710 cells
clusterName=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony
soClustering=RNA_snn_res.0.2

$Rscript $integration --Project $project --clusterName $clusterName --SeuratObject $seuratObject --SoClustering $soClustering --SoIdentities $soIdentities --ReducedDims $reducedDims --UMAP $umap --resolutions $resolutions --plotDir $plotDir --PredictionSuffix $PredictionSuffix --ptSize $ptSize --nCells_ATAC $nCells_ATAC --nCells_RNA $nCells_RNA

#### General statistics resulting from the integration
stats=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/integrationStats.R
project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k
resolutions=0.1,0.2,0.4,0.6,0.8
plotDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/integration_roundTwo/myeloid/
PredictionSuffix=Un_Round2
clusterName=myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony

$Rscript $stats --Project $project --resolutions $resolutions --plotDir $plotDir --suffix $PredictionSuffix --clusterName $clusterName

#### Add annotation column step ommited since we're using results of integration for rest of analyses
#### Removing MAST cells cluster

$Rscript /mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/scripts/removeMASTcluster_myeloid_Feb15_2023.R

#### Add peaks parent cell type level (pseudobulk is made using all cells)
addPeaks=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/addPeaks_byGroup.R
bedToBigBed=~/bedToBigBed
bedSort=~/bedSort
chromSizes=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/archR_aggr/peaks/pbs/hg38.chrom.sizes
groupColName=celltype
project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_noMAST
outPeakDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/peaks_Feb15_2023/myeloid_noMAST/

$Rscript $addPeaks --Project $project --groupColName $groupColName --outDir $outPeakDir --chromSizes $chromSizes --bedSort $bedSort --bedToBigBed $bedToBigBed

#### Add peaks subpopulation level (one pseudobulk per each subpopulation)
addPeaks=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/addPeaks_byGroup.R
bedToBigBed=~/bedToBigBed
bedSort=~/bedSort
chromSizes=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/archR_aggr/peaks/pbs/hg38.chrom.sizes

project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_noMAST
groupColName=predictedGroup_Un
outPeakDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/peaks_byCelltype/02_13_23/myeloid/
$Rscript $addPeaks --Project $project --groupColName $groupColName --outDir $outPeakDir --chromSizes $chromSizes --bedSort $bedSort --bedToBigBed $bedToBigBed

#### Make bigWig files parent cell type level (pseudobulk is made using all cells) and subpopulation level (one pseudobulk per each subpopulation)
script=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/archR/scripts/getGroupBw.R
Project=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k_noMAST
outDir=/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/UCSC_tracks_archR/bigwigs_Feb13_2023/myeloid
groupBy=celltype,predictedGroup_Un
tileSize=10
normMethod=ReadsInTSS,ReadsInPromoter,nFrags,nCells

$Rscript $script --Project $Project --outDir $outDir --groupBy $groupBy --normMethod $normMethod --tileSize $tileSize
