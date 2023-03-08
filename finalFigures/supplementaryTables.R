# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Merging ArchR summaries from projects using single and specific libraries in a single file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


paths <- c('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/B/B_R24DICE_Hu_B1_Bcell_ATAC_TSS8_nFragslog3.5_varFeats30k',
'/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/NK/NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k',
'/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/myeloid/myeloid_R24DICE_Hu_M2_mye_ATAC_TSS8_nFragslog3.5_varFeats30k',
'/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/lungCD4/lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k',
'/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/plots/lungCD8/lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')

paths <- paste0(paths, '/projectSummary.csv')

df <- paths %>% lapply(read.csv) %>% bind_rows
df$cellType <- c('B', 'NK', 'Myeloid', 'LungCD4', 'LungCD8')

write.csv(df, '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/tables4paper/ArchRStats.csv', row.names = FALSE, quote = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Writing cells per cluster tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lungCD4 <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')
lungCD8 <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k')
NK <- loadArchRProject('/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/projects/NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k')

lungCD4.table <- as.data.frame(table(lungCD4$lungCD4_R24_P_Hu_2_CD4_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.8))
names(lungCD4.table) <- c('Cluster', 'nCells')
lungCD4.table <- lungCD4.table %>% mutate(perc = nCells / sum(nCells))

lungCD8.table <- as.data.frame(table(lungCD8$lungCD8_R24_P_Hu_3_CD8_4D_ATAC_TSS10_nFragslog3.5_varFeats30k_harmony_0.4))
names(lungCD8.table) <- c('Cluster', 'nCells')
lungCD8.table <- lungCD8.table %>% mutate(perc = nCells / sum(nCells))

NK.table <- as.data.frame(table(NK$NK_R24DICE_Hu_N1_NK_ATAC_TSS8_nFragslog3.5_varFeats30k_harmony_0.4))
names(NK.table) <- c('Cluster', 'nCells')
NK.table <- NK.table %>% mutate(perc = nCells / sum(nCells))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Writing cells per subpopulation in each cell type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lungCD4.table.sp <- as.data.frame(table(lungCD4$annotations_12_23_22))
names(lungCD4.table.sp) <- c('subpopulation', 'nCells')
lungCD4.table.sp <- lungCD4.table.sp %>% mutate(perc = (nCells / sum(nCells)))
lungCD8.table.sp <- as.data.frame(table(lungCD8$annotations_12_23_22))
names(lungCD8.table.sp) <- c('subpopulation', 'nCells')
lungCD8.table.sp <- lungCD8.table.sp %>% mutate(perc = (nCells / sum(nCells)))
myeloid.table.sp <- as.data.frame(table(myeloid$predictedGroup_Un))
names(myeloid.table.sp) <- c('subpopulation', 'nCells')
myeloid.table.sp <- myeloid.table.sp %>% mutate(perc = (nCells / sum(nCells)))
NK.table.sp <- as.data.frame(table(NK$annotations_12_23_22))
names(NK.table.sp) <- c('subpopulation', 'nCells')
NK.table.sp <- NK.table.sp %>% mutate(perc = (nCells / sum(nCells)))

write.csv(lungCD4.table.sp, '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/tables4paper/cells_perSubpopulation_Mar1_2023/lungCD4.csv', row.names = FALSE, quote = FALSE)
write.csv(lungCD8.table.sp, '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/tables4paper/cells_perSubpopulation_Mar1_2023/lungCD8.csv', row.names = FALSE, quote = FALSE)
write.csv(myeloid.table.sp, '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/tables4paper/cells_perSubpopulation_Mar1_2023/myeloid.csv', row.names = FALSE, quote = FALSE)
write.csv(NK.table.sp, '/mnt/bioadhoc-temp/Groups/vd-vijay/acastillo/R24/singleLane_4paper/tables4paper/cells_perSubpopulation_Mar1_2023/NK.csv', row.names = FALSE, quote = FALSE)
