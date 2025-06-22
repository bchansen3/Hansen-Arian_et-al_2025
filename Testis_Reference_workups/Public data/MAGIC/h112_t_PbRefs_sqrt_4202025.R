
library(dplyr)
library(Seurat)
library(data.table)
library(tidyverse)


h_adult_112 <- readRDS("~/Documents/h_adult_112.rds")


a <- tibble(cellID = colnames(h_adult_112), clusterID = Idents(h_adult_112))

table(a$clusterID)

#a = object
#b = downsample size
#c = seed set
samp_ha <- function(a, b, c){
  
  lg <- subset(a, celltype == "Leydig")
  d_spg <- subset(a, celltype == "d_spg")
  endo <- subset(a, celltype == "Endo")
  ptm <- subset(a, celltype == "PTM")
  ssc <- subset(a, celltype == "SSCs")
  srt <- subset(a, celltype == "Sertoli")
  spc <- subset(a, celltype == "meiosis_spcs")
  imm <- subset(a, celltype == "Immune")
  
  set.seed(c)
  ds_lg <- if(table(Idents(lg)) > b){
            lg[, sample(colnames(lg), size = b, replace = F)]
  } else {
      lg
    }
  set.seed(c)
  ds_d_spg <- if(table(Idents(d_spg)) > b){
            d_spg[, sample(colnames(d_spg), size = b, replace = F)]
  } else {
      d_spg
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
  set.seed(c)
  ds_spc <- if(table(Idents(spc)) > b){
            spc[, sample(colnames(spc), size = b, replace = F)]
  } else {
      spc
    }
  set.seed(c)
  ds_imm <- if(table(Idents(imm)) > b){
            imm[, sample(colnames(imm), size = b, replace = F)]
  } else {
      imm
    }
  
  merger <- merge(x = ds_lg, y = list(ds_d_spg, ds_endo, ds_ptm, ds_ssc, ds_srt, ds_spc, ds_imm), merge.data = T)
  return(merger)
}

df_samp <- samp_ha(h_adult_112, 150, 23)

samp2 <- function(a, b, c, d){
  ct <- subset(a, idents = d)
  set.seed(c)
  subs <- ct[, sample(colnames(ct), size = b, replace = F)]
  return(subs)
}

#removing the sampled cells to make pseudobulk from the remainder
cells_for_pseudo <- colnames(h_adult_112)[!(colnames(h_adult_112) %in% colnames(df_samp))]
# re-add the Sertoli cells
sertoli_cells <- colnames(subset(h_adult_112, celltype =="Sertoli"))
cells_for_pseudo <- union(cells_for_pseudo, sertoli_cells)




magic_h112 <- fread("~/Documents/magic_sqrt_updated_41925/ha112_MAGIC_matrix_SQRT.csv")%>% t()
colnames(magic_h112) <- magic_h112[1,]
magic_h112 <- magic_h112[-1,]


mgc_pseudo <- as.matrix(magic_h112[,cells_for_pseudo])

pseudo <- h_adult_112[, cells_for_pseudo]
pseudo[["mgc"]] <- CreateAssayObject(data = mgc_pseudo)


table(pseudo$celltype)

#      Leydig      cnd_spm        d_spg        spmtd         Endo          PTM         SSCs meiosis_spcs       Immune 
#         479         2404          288          970          202          154          141          559          194 
#     Sertoli 
#          49 


#We will be making 5 pseudobulks total with 500 cells total - refer to code for ground truth values of cell type proportions



#1 90% Leydig cells
lg <- data.frame(samp2(pseudo, 450, 23, "Leydig")@assays$mgc$data)
dspg <- data.frame(samp2(pseudo, 10, 23, "d_spg")@assays$mgc$data)
ptm <- data.frame(samp2(pseudo, 10, 23, "PTM")@assays$mgc$data)
ssc <- data.frame(samp2(pseudo, 10, 23, "SSCs")@assays$mgc$data)
spc <- data.frame(samp2(pseudo, 10, 23, "meiosis_spcs")@assays$mgc$data)
srt <- data.frame(samp2(pseudo, 10, 23, "Sertoli")@assays$mgc$data)

pseudobulk <- cbind(lg,dspg,ptm,ssc,spc,srt) %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%  # Convert everything to numeric if necessary
  rowMeans(na.rm = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file = "~/Documents/MAGIC_T_04202025/h112_out/adult112_magic_lg90.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 75% Leydig cells
lg <- data.frame(samp2(pseudo, 375, 23, "Leydig")@assays$mgc$data)
dspg <- data.frame(samp2(pseudo, 25, 23, "d_spg")@assays$mgc$data)
ptm <- data.frame(samp2(pseudo, 25, 23, "PTM")@assays$mgc$data)
ssc <- data.frame(samp2(pseudo, 25, 23, "SSCs")@assays$mgc$data)
spc <- data.frame(samp2(pseudo, 25, 23, "meiosis_spcs")@assays$mgc$data)
srt <- data.frame(samp2(pseudo, 25, 23, "Sertoli")@assays$mgc$data)

pseudobulk <- cbind(lg,dspg,ptm,ssc,spc,srt) %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%  # Convert everything to numeric if necessary
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file = "~/Documents/MAGIC_T_04202025/h112_out/adult112_magic_lg75.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 50% Leydig cells
lg <- data.frame(samp2(pseudo, 250, 23, "Leydig")@assays$mgc$data)
dspg <- data.frame(samp2(pseudo, 50, 23, "d_spg")@assays$mgc$data)
ptm <- data.frame(samp2(pseudo, 50, 23, "PTM")@assays$mgc$data)
ssc <- data.frame(samp2(pseudo, 51, 23, "SSCs")@assays$mgc$data)
spc <- data.frame(samp2(pseudo, 50, 23, "meiosis_spcs")@assays$mgc$data)
srt <- data.frame(samp2(pseudo, 49, 23, "Sertoli")@assays$mgc$data)

pseudobulk <- cbind(lg,dspg,ptm,ssc,spc,srt) %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%  # Convert everything to numeric if necessary
  rowMeans(na.rm = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"



write.table(pseudobulk, file = "~/Documents/MAGIC_T_04202025/h112_out/adult112_magic_lg50.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 25% Leydig cells
lg <- data.frame(samp2(pseudo, 125, 23, "Leydig")@assays$mgc$data)
dspg <- data.frame(samp2(pseudo, 81, 23, "d_spg")@assays$mgc$data)
ptm <- data.frame(samp2(pseudo, 81, 23, "PTM")@assays$mgc$data)
ssc <- data.frame(samp2(pseudo, 82, 23, "SSCs")@assays$mgc$data)
spc <- data.frame(samp2(pseudo, 82, 23, "meiosis_spcs")@assays$mgc$data)
srt <- data.frame(samp2(pseudo, 49, 23, "Sertoli")@assays$mgc$data)

pseudobulk <- cbind(lg,dspg,ptm,ssc,spc,srt) %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%  # Convert everything to numeric if necessary
  rowMeans(na.rm = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file =  "~/Documents/MAGIC_T_04202025/h112_out/adult112_magic_lg25.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#5 10% Leydig cells
lg <- data.frame(samp2(pseudo, 50, 23, "Leydig")@assays$mgc$data)
dspg <- data.frame(samp2(pseudo, 101, 23, "d_spg")@assays$mgc$data)
ptm <- data.frame(samp2(pseudo, 100, 23, "PTM")@assays$mgc$data)
ssc <- data.frame(samp2(pseudo, 100, 23, "SSCs")@assays$mgc$data)
spc <- data.frame(samp2(pseudo, 100, 23, "meiosis_spcs")@assays$mgc$data)
srt <- data.frame(samp2(pseudo, 49, 23, "Sertoli")@assays$mgc$data)

pseudobulk <- cbind(lg,dspg,ptm,ssc,spc,srt) %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%  # Convert everything to numeric if necessary
  rowMeans(na.rm = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "~/Documents/MAGIC_T_04202025/h112_out/adult112_magic_lg10.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

### create the df_samp MAGIC object

ref_columns <- colnames(df_samp)

magic_ref <- magic_h112[,ref_columns]

df_samp[["mgc"]] <- CreateAssayObject(counts = magic_ref)

### MAKE CELL REFERENCE
a <- tibble(cellID = colnames(df_samp), clusterID = Idents(df_samp))
b <- data.frame(GetAssayData(df_samp, assay = "mgc"))
b$genes <- rownames(b)
#sum(duplicated(b$genes))


human_to_rat <- readRDS("~/Documents/Rat_to_Human.rds")
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
            file = "~/Documents/MAGIC_T_04202025/h112_out/h_a_112_MAGIC_Ref.txt", 
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = F)

