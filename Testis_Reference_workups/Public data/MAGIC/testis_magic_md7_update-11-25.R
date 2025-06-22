library(dplyr)
library(Seurat)
library(data.table)
library(tidyverse)


m_neonate_d7 <- readRDS("~/Documents/m_neonateD7.rds")


m_neonate_d7 <- SetIdent(m_neonate_d7, value= 'celltype')

a <- tibble(cellID = colnames(m_neonate_d7), clusterID = Idents(m_neonate_d7))



table(a$clusterID)


m_n_d7_labels <- m_neonate_d7


cluster.ids <- c("Sertoli", "Stroma/Leydig", "Pro-SG/SSC", "PTM", "Other", "diff_SPG", "Stroma/Leydig", "Immune", "Endo", "Lymph") 

names(cluster.ids) <- levels(m_n_d7_labels)

m_n_d7_labels <- RenameIdents(m_n_d7_labels, cluster.ids)
m_n_d7_labels$celltype <- Idents(m_n_d7_labels)

m_n_d7_labels <- SetIdent(m_n_d7_labels, value= 'celltype')

a <- tibble(cellID = colnames(m_n_d7_labels), clusterID = Idents(m_n_d7_labels))
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

df_samp <- samp_m(m_n_d7_labels, 150, 23)

samp2 <- function(a, b, c, d){
  ct <- subset(a, idents = d)
  set.seed(c)
  subs <- ct[, sample(colnames(ct), size = b, replace = F)]
  return(subs)
}

#removing the sampled cells to make pseudobulk from the remainder
cells_for_pseudo <- colnames(m_n_d7_labels)[!(colnames(m_n_d7_labels) %in% colnames(df_samp))]

pseudo <- m_n_d7_labels[, cells_for_pseudo]

### Load magic data

magic_mn7 <- fread("~/Documents/magic_sqrt_updated_41925/mn7_MAGIC_matrix_SQRT.csv")%>% t()
colnames(magic_mn7) <- magic_mn7[1,]
magic_mn7 <- magic_mn7[-1,]


mn7_mgc_pseudo <- magic_mn7[,cells_for_pseudo]


pseudo[["magic"]] <- CreateAssayObject(data = mn7_mgc_pseudo)


sampmagic <- magic_mn7[, !(colnames(magic_mn7) %in% cells_for_pseudo)]

df_samp[["magic"]] <- CreateAssayObject(data = sampmagic)


###############################

#1 90% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 900, 23, "Sertoli")@assays$magic$data)
lgSt <- data.frame(samp2(pseudo, 33, 23, "Stroma/Leydig")@assays$magic$data)
ptm <- data.frame(samp2(pseudo, 33, 23, "PTM")@assays$magic$data)
ssc <- data.frame(samp2(pseudo, 34, 23, "Pro-SG/SSC")@assays$magic$data)


pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file = "~/Documents/MAGIC_T_04202025/mn7_out/m_neo_d7_MAGIC_srt90.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 75% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 750, 23, "Sertoli")@assays$magic$data)
lgSt <- data.frame(samp2(pseudo, 83, 23, "Stroma/Leydig")@assays$magic$data)
ptm <- data.frame(samp2(pseudo, 83, 23, "PTM")@assays$magic$data)
ssc <- data.frame(samp2(pseudo, 84, 23, "Pro-SG/SSC")@assays$magic$data)

pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"



write.table(pseudobulk, file = "~/Documents/MAGIC_T_04202025/mn7_out/m_neo_d7_MAGIC_srt75.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)



#3 50% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 500, 23, "Sertoli")@assays$magic$data)
lgSt <- data.frame(samp2(pseudo, 166, 23, "Stroma/Leydig")@assays$magic$data)
ptm <- data.frame(samp2(pseudo, 166, 23, "PTM")@assays$magic$data)
ssc <- data.frame(samp2(pseudo, 168, 23, "Pro-SG/SSC")@assays$magic$data)

pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file = "~/Documents/MAGIC_T_04202025/mn7_out/m_neo_d7_MAGIC_srt50.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#4 25% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 250, 23, "Sertoli")@assays$magic$data)
lgSt <- data.frame(samp2(pseudo, 250, 23, "Stroma/Leydig")@assays$magic$data)
ptm <- data.frame(samp2(pseudo, 250, 23, "PTM")@assays$magic$data)
ssc <- data.frame(samp2(pseudo, 250, 23, "Pro-SG/SSC")@assays$magic$data)

pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file =  "~/Documents/MAGIC_T_04202025/mn7_out/m_neo_d7_MAGIC_srt25.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#5 10% Sertoli cells use 1k cells
srt <- data.frame(samp2(pseudo, 100, 23, "Sertoli")@assays$magic$data)
lgSt <- data.frame(samp2(pseudo, 300, 23, "Stroma/Leydig")@assays$magic$data)
ptm <- data.frame(samp2(pseudo, 300, 23, "PTM")@assays$magic$data)
ssc <- data.frame(samp2(pseudo, 300, 23, "Pro-SG/SSC")@assays$magic$data)

pseudobulk <- cbind(srt,lgSt,ptm,ssc) %>%
  mutate_all(as.numeric) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "~/Documents/MAGIC_T_04202025/mn7_out/m_neo_d7_MAGIC_srt10.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

### MAKE CELL REFERENCE
a <- tibble(cellID = colnames(df_samp), clusterID = Idents(df_samp))
b <- data.frame(GetAssayData(df_samp, assay = "magic"))
b$genes <- rownames(b)
#sum(duplicated(b$genes))


# Keep only the Gene Names
mouse_to_rat <-  readRDS("~/Downloads/Rat_to_Mouse.rds")

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
            file = "~/Documents/MAGIC_T_04202025/mn7_out/m_neo_d7_MAGIC.txt", 
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = F)

