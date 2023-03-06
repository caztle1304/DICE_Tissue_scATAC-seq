# scATAC-seq DICE Tissue
--------------
- Angel Adrian De la Cruz Castillo (acastillo@lji.org)
- Vijayanand Lab (https://www.lji.org/labs/vijayanand/)
- La Jolla Institute for Immunology (LJI)
- La Jolla, CA, USA

## Description
--------------
This repository contains scripts used to analyze scATAC-seq data ford DICE Tissue project (unpublished). It contains R scripts that were called in command line

## Requirements 
--------------
- R (v4.0.1)
- ArchR (v1.0.1)
- Macs2 (v2.2.7.1)
- Cell Ranger ATAC (v2.0.0)
- bedToBigBed (v2.8)
- Chrom sizes file used in bigBed file formatting was obtained using fetchChromSizes from UCSC utilities with the following command line:

```console
fetchChromSizes=~/fetchChromSizes
$fetchChromSizes hg38 hg38.chrom.sizes
```
- bedSort from UCSC utilities was used as well to format bigBed files

To download UCSC utilities go to: https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

## Raw data pre-processing 
--------------
- Pre-processing was done entirely using Cell Ranger ATAC mkffastq and count (v2.0.0)
- sh files contain lines and parameters that were used to call R scripts to process fragment files and visualize results. First run_createProjects.sh was used for creation of all initial ArchR projects and then, the rest of sh files were used for each cell-type-specific ArchR project. 
- Scripts used to create figures in paper are located in finalFigures folder

## Contact 
--------------
Please email Angel De la Cruz Castillo (acastillo@lji.org) and/or Pandurangan Vijayanand(vijay@lji.org)
