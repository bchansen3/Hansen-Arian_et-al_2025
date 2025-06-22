
library(dplyr)
library(Seurat)
library(data.table)
library(tidyverse)


h_infant_124 <- readRDS("~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/objectsforcounts_50825/h_infant_124.rds")


h_infant_124 <- SetIdent(h_infant_124, value= "celltype")
new_names <- c("Sertoli", "Leydig", "PTM", "Germ", "Endo")
names(new_names) <- levels(h_infant_124)
h_infant_124 <- RenameIdents(object = h_infant_124, new_names)
h_infant_124$celltype <- h_infant_124@active.ident



a <- tibble(cellID = colnames(h_infant_124), clusterID = Idents(h_infant_124))
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
  
  
  merger <- merge(x = ds_lg, y = list(ds_germ, ds_endo, ds_ptm, ds_srt), merge.data = T)
  return(merger)
}

df_samp <- samp_ha(h_infant_124, 100, 23)

df_samp <- JoinLayers(df_samp)


samp2 <- function(a, b, c, d){
  ct <- subset(a, idents = d)
  set.seed(c)
  subs <- ct[, sample(colnames(ct), size = b, replace = F)]
  return(subs)
}

#removing the sampled cells to make pseudobulk from the remainder
cells_for_pseudo <- colnames(h_infant_124)[!(colnames(h_infant_124) %in% colnames(df_samp))]

pseudo <- h_infant_124[, cells_for_pseudo]

table(pseudo$celltype)



#We will be making 5 pseudobulks total with 500 cells total - refer to code for ground truth values of cell type proportions



#1 90% Sertoli cells (250 cells)
lg <- data.frame(samp2(pseudo, 12, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 13, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 13, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 12, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 450, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  mutate(across(everything(), as.numeric)) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h124/H_I-124_COUNTS_SRT90.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 75% Sertoli cells (252 cells)
lg <- data.frame(samp2(pseudo, 50, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 50, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 50, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 50, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 450, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  mutate(across(everything(), as.numeric)) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h124/H_I-124_COUNTS_SRT75.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 50% Sertoli cells
lg <- data.frame(samp2(pseudo, 63, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 62, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 63, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 62, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 250, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  mutate(across(everything(), as.numeric)) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h124/H_I-124_COUNTS_SRT50.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 25% Sertoli cells 252 cells
lg <- data.frame(samp2(pseudo, 99, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 80, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 98, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 98, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 125, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  mutate(across(everything(), as.numeric)) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h124/H_I-124_COUNTS_SRT25.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

#5 10% Sertoli cells
lg <- data.frame(samp2(pseudo, 124, 23, "Leydig")@assays$RNA$counts)
germ <- data.frame(samp2(pseudo, 80, 23, "Germ")@assays$RNA$counts)
ptm <- data.frame(samp2(pseudo, 123, 23, "PTM")@assays$RNA$counts)
endo <- data.frame(samp2(pseudo, 123, 23, "Endo")@assays$RNA$counts)
srt <- data.frame(samp2(pseudo, 50, 23, "Sertoli")@assays$RNA$counts)

pseudobulk <- cbind(lg,germ,ptm,endo,srt) %>%
  mutate(across(everything(), as.numeric)) %>%
  rowSums() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"


write.table(pseudobulk, file =  "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/PBs/h124/H_I-124_COUNTS_SRT10.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)



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
colnames(b) <- gsub("\\.", "-", colnames(b))


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
            file = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/raw_counts/TESTIS/h_i_124_MAGIC_REF.txt", 
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = F)

