## Cell ranger lines used to demultiplex and generate fragments file

#### Lung CD4 #####

# mkfastq
/mnt/BioAdHoc/Groups/vd-vijay/jrocha/cellranger/atac-seq/cellranger_atac_v2/cellranger-atac-2.0.0/bin/atac/mkfastq --run=/home/jrocha/BioAdHoc/scATAC-seq/Projects/R24/NV047/snake/210503_A00475_0311_BHKT23DSXY_NV047 --samplesheet=/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/CD4_CD8_cellrangerV2/NV047/NV047/data/Sample-sheet.csv --id=mkfastq --jobmode=torque --use-bases-mask=Y*,I8n*,Y16,Y*

# count
/home/jrocha/BioAdHoc/cellranger/atac-seq/cellranger_atac_v2/cellranger-atac-2.0.0/cellranger-atac count --id R24P_Hu_2_CD4_4D_ATAC --reference /home/jrocha/BioAdHoc/cellranger/atac-seq/cellranger_atac_v2/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --fastqs /mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/CD4_CD8_cellrangerV2/NV047/NV047/mkfastq/outs/fastq_path/R24_P_ATAC/R24P_Hu_2_CD4_4D_ATAC --jobmode torque --maxjobs 100

#### Lung CD8 ####

# mkfastq
/mnt/BioAdHoc/Groups/vd-vijay/jrocha/cellranger/atac-seq/cellranger_atac_v2/cellranger-atac-2.0.0/bin/atac/mkfastq --run=/home/jrocha/BioAdHoc/scATAC-seq/Projects/R24/NV047/snake/210503_A00475_0311_BHKT23DSXY_NV047 --samplesheet=/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/CD4_CD8_cellrangerV2/NV047/NV047/data/Sample-sheet.csv --id=mkfastq --jobmode=torque --use-bases-mask=Y*,I8n*,Y16,Y*

# count
/home/jrocha/BioAdHoc/cellranger/atac-seq/cellranger_atac_v2/cellranger-atac-2.0.0/cellranger-atac count --id R24P_Hu_3_CD8_4D_ATAC --reference /home/jrocha/BioAdHoc/cellranger/atac-seq/cellranger_atac_v2/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --fastqs /mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/CD4_CD8_cellrangerV2/NV047/NV047/mkfastq/outs/fastq_path/R24_P_ATAC/R24P_Hu_3_CD8_4D_ATAC --jobmode torque --maxjobs 100

#### Myeloid ####
# mkfastq
/mnt/BioAdHoc/Groups/vd-vijay/jrocha/cellranger/atac-seq/cellranger_atac_v2/cellranger-atac-2.0.0/bin/atac/mkfastq --run=/mnt/NovaSeq/211108_A00475_0377_AHTGK3DSX2_NV068 --samplesheet=/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/NV068/NV068_cellRangerV2/NV068/data/Sample-sheet.csv --id=mkfastq --jobmode=torque --use-bases-mask=Y*,I8n*,Y16,Y*

# count
/home/jrocha/BioAdHoc/cellranger/atac-seq/cellranger_atac_v2/cellranger-atac-2.0.0/cellranger-atac count --id=R24DICE_Hu_M2_mye_ATAC --fastqs=/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/NV068/NV068_cellRangerV2/NV068/mkfastq/outs/fastq_path/R24DICE_ATAC/R24DICE_Hu_M2_mye_ATAC --reference=/home/jrocha/BioAdHoc/cellranger/atac-seq/cellranger_atac_v2/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --jobmode=torque

#### NK ####
# mkfastq
/mnt/BioAdHoc/Groups/vd-vijay/jrocha/cellranger/atac-seq/cellranger_atac_v2/cellranger-atac-2.0.0/bin/atac/mkfastq --run=/mnt/NovaSeq/211108_A00475_0377_AHTGK3DSX2_NV068 --samplesheet=/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/NV068/NV068_cellRangerV2/NV068/data/Sample-sheet.csv --id=mkfastq --jobmode=torque --use-bases-mask=Y*,I8n*,Y16,Y*

# count
/home/jrocha/BioAdHoc/cellranger/atac-seq/cellranger_atac_v2/cellranger-atac-2.0.0/cellranger-atac count --id=R24DICE_Hu_N1_NK_ATAC --fastqs=/mnt/BioAdHoc/Groups/vd-vijay/acastillo/cellRanger/NV068/NV068_cellRangerV2/NV068/mkfastq/outs/fastq_path/R24DICE_ATAC/R24DICE_Hu_N1_NK_ATAC --reference=/home/jrocha/BioAdHoc/cellranger/atac-seq/cellranger_atac_v2/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --jobmode=torque
