
library(dplyr)
library(Seurat)
library(data.table)


m_neonate_d2 <- readRDS("~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/objectsforcounts_50825/m_neonateD2.rds")


m_neonate_d2 <- SetIdent(m_neonate_d2, value= 'celltype')

a <- tibble(cellID = colnames(m_neonate_d2), clusterID = Idents(m_neonate_d2))



table(a$clusterID)

#    Sertoli     Stroma/Leydig     PTM        Pro-SG/SSC        Endo          Immune 
#     8887          6022          2065          1329           219           164 

# Just use the Sertoli     Stroma/Leydig     PTM        Pro-SG/SSC

#a = object
#b = downsample size
#c = seed set
samp_m <- function(a, b, c){
  
  lgSt <- subset(a, celltype == "Stroma/Leydig")
  endo <- subset(a, celltype == "Endo")
  ptm <- subset(a, celltype == "PTM")
  srt <- subset(a, celltype == "Sertoli")
  ssc <- subset(a, celltype == "Pro-SG/SSC")
  imm <- subset(a, celltype == "Immune")
  
  set.seed(c)
  ds_lgSt <- if(table(Idents(lgSt)) > b){
    lgSt[, sample(colnames(lgSt), size = b, replace = F)]
  } else {
    lgSt
  }
  set.seed(c)
  ds_ptm <- if(table(Idents(ptm)) > b){
    ptm[, sample(colnames(ptm), size = b, replace = F)]
  } else {
    ptm
  }
  set.seed(c)
  ds_ssc <- if(table(Idents(ssc)) > b){
    ssc[, sample(colnames(ssc), size = b, replace = F)]
  } else {
    ssc
  }
  set.seed(c)
  ds_srt <- if(table(Idents(srt)) > b){
    srt[, sample(colnames(srt), size = b, replace = F)]
  } else {
    srt
  }
  merger <- merge(x = ds_lgSt, y = list( ds_ptm, ds_ssc, ds_srt), merge.data = T)
  return(merger)
}

df_samp <- samp_m(m_neonate_d2, 150, 23)
df_samp <- JoinLayers(df_samp)

samp2 <- function(a, b, c, d){
  ct <- subset(a, idents = d)
  set.seed(c)
  subs <- ct[, sample(colnames(ct), size = b, replace = F)]
  return(subs)
}

#removing the sampled cells to make pseudobulk from the remainder
cells_for_pseudo <- colnames(m_neonate_d2)[!(colnames(m_neonate_d2) %in% colnames(df_samp))]

pseudo <- m_neonate_d2[, cells_for_pseudo]

table(pseudo$celltype)


###############################

#1 90% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 900, 23, "Sertoli")@assays$RNA$counts)
lgSt <- data.frame(samp2(pseudo, 33, 23, "Stroma/Leydig")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 33, 23, "PTM")@assays$RNA$counts)
ssc <- data.frame(samp2(pseudo, 34, 23, "Pro-SG/SSC")@assays$RNA$counts)


pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/md2/m_neo_d2_COUNTS_srt90.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 75% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 750, 23, "Sertoli")@assays$RNA$counts)
lgSt <- data.frame(samp2(pseudo, 83, 23, "Stroma/Leydig")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 83, 23, "PTM")@assays$RNA$counts)
ssc <- data.frame(samp2(pseudo, 84, 23, "Pro-SG/SSC")@assays$RNA$counts)

pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"



write.table(pseudobulk, file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/md2/m_neo_d2_COUNTS_srt75.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)



#3 50% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 500, 23, "Sertoli")@assays$RNA$counts)
lgSt <- data.frame(samp2(pseudo, 166, 23, "Stroma/Leydig")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 166, 23, "PTM")@assays$RNA$counts)
ssc <- data.frame(samp2(pseudo, 168, 23, "Pro-SG/SSC")@assays$RNA$counts)

pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/md2/m_neo_d2_COUNTS_srt50.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#4 25% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 250, 23, "Sertoli")@assays$RNA$counts)
lgSt <- data.frame(samp2(pseudo, 250, 23, "Stroma/Leydig")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 250, 23, "PTM")@assays$RNA$counts)
ssc <- data.frame(samp2(pseudo, 250, 23, "Pro-SG/SSC")@assays$RNA$counts)

pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/md2/m_neo_d2_COUNTS_srt25.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#5 10% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 100, 23, "Sertoli")@assays$RNA$counts)
lgSt <- data.frame(samp2(pseudo, 300, 23, "Stroma/Leydig")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 300, 23, "PTM")@assays$RNA$counts)
ssc <- data.frame(samp2(pseudo, 300, 23, "Pro-SG/SSC")@assays$RNA$counts)

pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/md2/m_neo_d2_COUNTS_srt10.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

### MAKE CELL REFERENCE
a <- tibble(cellID = colnames(df_samp), clusterID = Idents(df_samp))
b <- data.frame(GetAssayData(df_samp, assay = "RNA", layer = "counts"))
b$genes <- rownames(b)
#sum(duplicated(b$genes))


# Keep only the Gene Names
mouse_to_rat <-  readRDS("~/OneDrive - UW/data for deconv/public data (B)/Testis/Rat_to_Mouse.rds")

mouse_to_rat <- mouse_to_rat[,c("Mouse_GRCm39_GeneSymbol","Rat_mRatBN7.2_GeneSymbol")]
colnames(mouse_to_rat) <- c("mouse","rat")
mouse_to_rat <- mouse_to_rat %>% distinct() 
mouse_to_rat <- mouse_to_rat %>%
  add_count(rat) %>%           # Add count of occurrences for 'rat'
  filter(n == 1) %>%           # Keep only rows where the 'rat' count is 1 (no duplicates)
  select(-n) 

b <- left_join(b,mouse_to_rat, by = c('genes' = 'mouse')) %>% group_by(rat) %>%
  filter(n() == 1) %>%  ungroup() %>% column_to_rownames(., var="rat") %>% select(-genes)


names(b) <- gsub("\\.", "-", names(b))


t_a <- data.frame(t(a))
colnames(t_a) <- t_a[1,]

cell_ref <- rbind(t_a, b)
table(colnames(cell_ref) == cell_ref["cellID",])

cell_ref$GENEID <- rownames(cell_ref)
cell_ref <- cell_ref %>%
  relocate(GENEID)

#Getting rid of the cellID row so the cluster ID is now at the top
cell_ref <- cell_ref[-1,]

write.table(cell_ref, 
            file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/m_neo_d2_COUNTS_CELLREF.txt", 
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = F)


