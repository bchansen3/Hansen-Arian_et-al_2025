
library(dplyr)
library(Seurat)
library(data.table)
library(tidyverse)


h_infant_120 <- readRDS("~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/objectsforcounts_50825/h_infant_120.rds")

h_infant_120 <- SetIdent(h_infant_120, value= 'celltype')

a <- tibble(cellID = colnames(h_infant_120), clusterID = Idents(h_infant_120))

table(a$clusterID)


#a = object
#b = downsample size
#c = seed set
samp_ha <- function(a, b, c){
  
  lg <- subset(a, celltype == "Leydig")
  germ <- subset(a, celltype == "Germ")
  endo <- subset(a, celltype == "Endo")
  ptm <- subset(a, celltype == "PTM")
  srt <- subset(a, celltype == "Sertoli")
  imm <- subset(a, celltype == "Immune")
  
  set.seed(c)
  ds_lg <- if(table(Idents(lg)) > b){
            lg[, sample(colnames(lg), size = b, replace = F)]
  } else {
      lg
    }
  set.seed(c)
  ds_germ <- if(table(Idents(germ)) > b){
            germ[, sample(colnames(germ), size = b, replace = F)]
  } else {
      germ
    }
  set.seed(c)
  ds_endo <- if(table(Idents(endo)) > b){
            endo[, sample(colnames(endo), size = b, replace = F)]
  } else {
      endo
    }
  set.seed(c)
  ds_ptm <- if(table(Idents(ptm)) > b){
            ptm[, sample(colnames(ptm), size = b, replace = F)]
  } else {
      ptm
    }
  set.seed(c)
  ds_srt <- if(table(Idents(srt)) > b){
            srt[, sample(colnames(srt), size = b, replace = F)]
  } else {
      srt
    }
  set.seed(c)
  ds_imm <- if(table(Idents(imm)) > b){
            imm[, sample(colnames(imm), size = b, replace = F)]
  } else {
      imm
    }
  
  merger <- merge(x = ds_lg, y = list(ds_germ, ds_endo, ds_ptm, ds_srt, ds_imm), merge.data = T)
  return(merger)
}

df_samp <- samp_ha(h_infant_120, 50, 23)

df_samp <- JoinLayers(df_samp)


samp2 <- function(a, b, c, d){
  ct <- subset(a, idents = d)
  set.seed(c)
  subs <- ct[, sample(colnames(ct), size = b, replace = F)]
  return(subs)
}

#removing the sampled cells to make pseudobulk from the remainder
cells_for_pseudo <- colnames(h_infant_120)[!(colnames(h_infant_120) %in% colnames(df_samp))]
# re-add the Sertoli cells
germ_cells <- colnames(subset(h_infant_120, celltype =="Germ"))
cells_for_pseudo <- union(cells_for_pseudo, germ_cells)


pseudo <- h_infant_120[, cells_for_pseudo]



#table(pseudo$celltype)

#      Leydig      cnd_spm        d_spg        spmtd         Endo          PTM         SSCs meiosis_spcs       Immune 
#         479         2404          288          970          202          154          141          559          194 
#     Sertoli 
#          49 


#We will be making 5 pseudobulks total with 500 cells total - refer to code for ground truth values of cell type proportions



#1 90% Sertoli cells (250 cells)
lg <- data.frame(samp2(pseudo, 6, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 7, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 6, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 6, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 225, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h120/H_I-120_COUNTS_SRT90.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 75% Sertoli cells (252 cells)
lg <- data.frame(samp2(pseudo, 16, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 16, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 16, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 15, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 189, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h120/H_I-120_COUNTS_SRT75.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 50% Sertoli cells
lg <- data.frame(samp2(pseudo, 31, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 32, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 31, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 31, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 125, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h120/H_I-120_COUNTS_SRT50.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 25% Sertoli cells 252 cells
lg <- data.frame(samp2(pseudo, 51, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 36, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 51, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 51, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 63, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h120/H_I-120_COUNTS_SRT25.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

#5 10% Sertoli cells
lg <- data.frame(samp2(pseudo, 63, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 36, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 63, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 63, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 25, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h120/H_I-120_COUNTS_SRT10.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)



### MAKE CELL REFERENCE
a <- tibble(cellID = colnames(df_samp), clusterID = Idents(df_samp))
b <- data.frame(GetAssayData(df_samp, assay = "RNA", layer = "counts"))
b$genes <- rownames(b)
#sum(duplicated(b$genes))


human_to_rat <- readRDS("~/OneDrive - UW/data for deconv/public data (B)/Testis/Rat_to_Human.rds")
# Keep only the Gene Names
human_to_rat <- human_to_rat[,c("Human_GRCh38_GeneSymbol","Rat_mRatBN7.2_GeneSymbol")]
colnames(human_to_rat) <- c("human","rat")
human_to_rat <- human_to_rat %>% distinct()

b <- left_join(b,human_to_rat, by = c('genes' = 'human')) %>% column_to_rownames(., var="rat") %>% select(-genes)


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
            file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h120/h_i_120_COUNTS_Ref.txt", 
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = F)

saveRDS(a,'~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/h120_testis_CellIds.rds')
