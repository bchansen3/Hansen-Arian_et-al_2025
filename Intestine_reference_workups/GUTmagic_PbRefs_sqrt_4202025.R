set.seed(23)


library(Seurat)
library(readr)
library(biomaRt)
library(tidyverse)
library(data.table)

#a = object
#b = downsample size
#c = seed set
samp <- function(a, b, c){
  
  ent <- subset(a, cell_type == "Enterocyte")
  gob <- subset(a, cell_type == "Goblet cell")
  other <- subset(a, cell_type == "Other")
  prol <- subset(a, cell_type == "Proliferative cell")
  best4 <- subset(a, cell_type == "BEST4")
  tuft <- subset(a, cell_type == "Tuft cell")
  eec <- subset(a, cell_type == "EEC")
  pan <- subset(a, cell_type == "Paneth cell")
  
  set.seed(c)
  dse <- if(table(Idents(ent)) > b){
    ent[, sample(colnames(ent), size = b, replace = F)]
  } else {
    ent
  }
  set.seed(c)
  dsg <- if(table(Idents(gob)) > b){
    gob[, sample(colnames(gob), size = b, replace = F)]
  } else {
    gob
  }
  set.seed(c)
  dso <- if(table(Idents(other)) > b){
    other[, sample(colnames(other), size = b, replace = F)]
  } else {
    other
  }
  set.seed(c)
  dsp <- if(table(Idents(prol)) > b){
    prol[, sample(colnames(prol), size = b, replace = F)]
  } else {
    prol
  }
  set.seed(c)
  dsb <- if(table(Idents(best4)) > b){
    best4[, sample(colnames(best4), size = b, replace = F)]
  } else {
    best4
  }
  set.seed(c)
  dst <- if(table(Idents(tuft)) > b){
    tuft[, sample(colnames(tuft), size = b, replace = F)]
  } else {
    tuft
  }
  set.seed(c)
  dsee <- if(table(Idents(eec)) > b){
    eec[, sample(colnames(eec), size = b, replace = F)]
  } else {
    eec
  }
  set.seed(c)
  dspan <- if(table(Idents(pan)) > b){
    pan[, sample(colnames(pan), size = b, replace = F)]
  } else {
    pan
  }
  
  merger <- merge(x = dse, y = list(dsg, dso, dsp, dsb, dst, dsee, dspan), merge.data = T)
  return(merger)
}

samp2 <- function(a, b, c, d){
  ct <- subset(a, idents = d)
  set.seed(c)
  subs <- ct[, sample(colnames(ct), size = b, replace = F)]
  return(subs)
}

samp_NOOTHER <- function(a, b, c){
  
  ent <- subset(a, cell_type == "Enterocyte")
  gob <- subset(a, cell_type == "Goblet cell")
  prol <- subset(a, cell_type == "Proliferative cell")
  best4 <- subset(a, cell_type == "BEST4")
  tuft <- subset(a, cell_type == "Tuft cell")
  eec <- subset(a, cell_type == "EEC")
  pan <- subset(a, cell_type == "Paneth cell")
  
  set.seed(c)
  dse <- if(table(Idents(ent)) > b){
    ent[, sample(colnames(ent), size = b, replace = F)]
  } else {
    ent
  }
  set.seed(c)
  dsg <- if(table(Idents(gob)) > b){
    gob[, sample(colnames(gob), size = b, replace = F)]
  } else {
    gob
  }
  set.seed(c)
  dsp <- if(table(Idents(prol)) > b){
    prol[, sample(colnames(prol), size = b, replace = F)]
  } else {
    prol
  }
  set.seed(c)
  dsb <- if(table(Idents(best4)) > b){
    best4[, sample(colnames(best4), size = b, replace = F)]
  } else {
    best4
  }
  set.seed(c)
  dst <- if(table(Idents(tuft)) > b){
    tuft[, sample(colnames(tuft), size = b, replace = F)]
  } else {
    tuft
  }
  set.seed(c)
  dsee <- if(table(Idents(eec)) > b){
    eec[, sample(colnames(eec), size = b, replace = F)]
  } else {
    eec
  }
  set.seed(c)
  dspan <- if(table(Idents(pan)) > b){
    pan[, sample(colnames(pan), size = b, replace = F)]
  } else {
    pan
  }
  
  merger <- merge(x = dse, y = list(dsg, dsp, dsb, dst, dsee, dspan), merge.data = T)
  return(merger)
}


setwd("~/Documents/magic_corrections")



#load in data
df <- readRDS("DS_paper1_cpm.rds")


#this is a matrix with cells as columns and rows as genes - cell names maintained. 
magic_df1 <- fread("FINAL_df1_MAGIC_matrix.csv")%>% t()
colnames(magic_df1) <- magic_df1[1,]
magic_df1 <- magic_df1[-1,]


magic_cells <- colnames(magic_df1)
df <- df[, magic_cells]
df[["MAGIC"]] <- CreateAssayObject(data = magic_df1)


#Need to make a "cell_type" column in metadata for below function to work
df$cell_type <- df@active.ident

#sample from data layer
df_samp <- samp_NOOTHER(df, 300, 23)

#removing the sampled cells to make pseudobulk from the remainder
cells_for_pseudo <- colnames(df)[!(colnames(df) %in% colnames(df_samp))]
pseudo <- df[, cells_for_pseudo]


"Proliferative cell         Enterocyte        Paneth cell 
               793               2892                541 "

#I'll do 1000 cells in each pseudobulk

#1 50% ents

ent <- data.frame(samp2(pseudo, 500, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 300, 23, "Proliferative cell")@assays$MAGIC$data)
pan <- data.frame(samp2(pseudo, 200, 23, "Paneth cell")@assays$MAGIC$data)


pseudobulk <- cbind(ent, prol, pan) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df1_MAGIC_pCb1.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2 10% ents ## doing 500 cells total
ent <- data.frame(samp2(pseudo, 100, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 700, 23, "Proliferative cell")@assays$MAGIC$data)
pan <- data.frame(samp2(pseudo, 200, 23, "Paneth cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, prol, pan) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df1_MAGIC_pCb2.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#3 25% ents ## 500 cells
ent <- data.frame(samp2(pseudo, 250, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 550, 23, "Proliferative cell")@assays$MAGIC$data)
pan <- data.frame(samp2(pseudo, 200, 23, "Paneth cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, prol, pan) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df1_MAGIC_pCb3.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#4 75% ents
ent <- data.frame(samp2(pseudo, 750, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 150, 23, "Proliferative cell")@assays$MAGIC$data)
pan <- data.frame(samp2(pseudo, 100, 23, "Paneth cell")@assays$MAGIC$data)


pseudobulk <- cbind(ent, prol, pan) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df1_MAGIC_pCb4.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#5 90% ents
ent <- data.frame(samp2(pseudo, 900, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 50, 23, "Proliferative cell")@assays$MAGIC$data)
pan <- data.frame(samp2(pseudo, 50, 23, "Paneth cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, prol, pan) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df1_MAGIC_pCb5.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


a <- tibble(cellID = colnames(df_samp), clusterID = Idents(df_samp))
b <- data.frame(df_samp@assays$MAGIC$data)
b <- b %>%  mutate_all(function(x) as.numeric((x)))

#b[c(1:10), c(1:10)]

colnames(b) <-gsub(".", "-", colnames(b), fixed = TRUE)
t_a <- data.frame(t(a))
colnames(t_a) <- t_a[1,]

cell_ref <- rbind(t_a, b)
#head(cell_ref)
table(colnames(cell_ref) == cell_ref["cellID",])

cell_ref$GENEID <- rownames(cell_ref)
cell_ref <- cell_ref %>%
  relocate(GENEID)

#Getting rid of the cellID row so the cluster ID is now at the top
cell_ref <- cell_ref[-1,]

write.table(cell_ref, 
            file = "Paper1_MAGIC_CellRef_final.txt", 
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)




##### PAPER 2


#load in data
df <- readRDS("DS_paper2_cpm.rds")

#load in saver data
#this is a matrix with cells as columns and rows as genes - cell names maintained. I'm going to remake the seurat so that the correct cells are present and add in the daver data from there
magic_df2 <- fread("FINAL_df2_MAGIC_matrix.csv")%>% t()
colnames(magic_df2) <- magic_df2[1,]
magic_df2 <- magic_df2[-1,]


magic_cells <- colnames(magic_df2)
df <- df[, magic_cells]
df[["MAGIC"]] <- CreateAssayObject(data = magic_df2)




#Need to make a "cell_type" column in metadata for below function to work
df$cell_type <- df@active.ident

#sample from data layer
df_samp <- samp(df, 300, 23)

#removing the sampled cells to make pseudobulk from the remainder
cells_for_pseudo <- colnames(df)[!(colnames(df) %in% colnames(df_samp))]
pseudo <- df[, cells_for_pseudo]



#1
ent <- data.frame(samp2(pseudo, 500, 23, "Enterocyte")@assays$MAGIC$data)
oth <- data.frame(samp2(pseudo, 200, 23, "Other")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 300, 23, "Goblet cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, oth, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, 
            file = "df2_MAGIC_pCb1.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2
ent <- data.frame(samp2(pseudo, 100, 23, "Enterocyte")@assays$MAGIC$data)
oth <- data.frame(samp2(pseudo, 400, 23, "Other")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 500, 23, "Goblet cell")@assays$MAGIC$data)


pseudobulk <- cbind(ent, oth, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df2_MAGIC_pCb2.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#3
ent <- data.frame(samp2(pseudo, 250, 23, "Enterocyte")@assays$MAGIC$data)
oth <- data.frame(samp2(pseudo, 250, 23, "Other")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 500, 23, "Goblet cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, oth, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df2_MAGIC_pCb3.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#4
ent <- data.frame(samp2(pseudo, 750, 23, "Enterocyte")@assays$MAGIC$data)
oth <- data.frame(samp2(pseudo, 100, 23, "Other")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 150, 23, "Goblet cell")@assays$MAGIC$data)


pseudobulk <- cbind(ent, oth, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df2_MAGIC_pCb4.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#5
ent <- data.frame(samp2(pseudo, 900, 23, "Enterocyte")@assays$MAGIC$data)
oth <- data.frame(samp2(pseudo, 50, 23, "Other")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 50, 23, "Goblet cell")@assays$MAGIC$data)


pseudobulk <- cbind(ent, oth, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

write.table(pseudobulk, file = "df2_MAGIC_pCb5.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


a <- tibble(cellID = colnames(df_samp), clusterID = Idents(df_samp))
b <- data.frame(df_samp@assays$MAGIC$data)
b <- b %>%  mutate_all(function(x) as.numeric((x)))



colnames(b) <-gsub(".", "-", colnames(b), fixed = TRUE)
t_a <- data.frame(t(a))
colnames(t_a) <- t_a[1,]

cell_ref <- rbind(t_a, b)
#head(cell_ref)
table(colnames(cell_ref) == cell_ref["cellID",])

cell_ref$GENEID <- rownames(cell_ref)
cell_ref <- cell_ref %>%
  relocate(GENEID)

#Getting rid of the cellID row so the cluster ID is now at the top
cell_ref <- cell_ref[-1,]

write.table(cell_ref, 
            file = "Paper2_MAGIC_CellRef_final.txt", 
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


## PAPER 3

#load in data
df <- readRDS("DS_paper3_cpm.rds")

#load in saver data
#this is a matrix with cells as columns and rows as genes - cell names maintained. I'm going to remake the seurat so that the correct cells are present and add in the daver data from there
magic_df3 <- fread("FINAL_df3_MAGIC_matrix.csv")%>% t()
colnames(magic_df3) <- magic_df3[1,]
magic_df3 <- magic_df3[-1,]


magic_cells <- colnames(magic_df3)
df <- df[, magic_cells]
df[["MAGIC"]] <- CreateAssayObject(data = magic_df3)


#Need to make a "cell_type" column in metadata for below function to work
df$cell_type <- df@active.ident

#sample from data layer
df_samp <- samp(df, 300, 23)
table(Idents(df_samp))

#removing the sampled cells to make pseudobulk from the remainder
cells_for_pseudo <- colnames(df)[!(colnames(df) %in% colnames(df_samp))]
pseudo <- df[, cells_for_pseudo]

#memory saver
rm(df)
gc()

library(biomaRt)

### RUN THIS FOR REDOS
#pseudo <- readRDS("~/OneDrive - UW/data for deconv/MAGIC_objects_CPM_Correction/pseudo_3_magic.rds")

# Proliferative cell         Enterocyte        Goblet cell              BEST4 
#           1565               1860               1006                 16 
#Lost a lot of cell types - I'll do 1000 cells in each pseudobulk and omit BEST4 cells because the count is so low

#I'll do 1000 cells in each pseudobulk

#1
ent <- data.frame(samp2(pseudo, 500, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 300, 23, "Proliferative cell")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 200, 23, "Goblet cell")@assays$MAGIC$data)


pseudobulk <- cbind(ent, prol, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

#need to rename genes from ensembl to gene names
ens <- pseudobulk$GENEID

mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
ens_to_symbol_biomart <- getBM(
  filters = 'ensembl_gene_id',
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  values = ens,
  mart = mart)

colnames(ens_to_symbol_biomart) <- c("GENEID", "name")

pseudobulk <- left_join(pseudobulk,ens_to_symbol_biomart, by="GENEID")

#removing any NAs, averaging duplicated genes
#for some reason there is a blank gene named "", I will remove this in the last line of the below code
pseudobulk <- pseudobulk %>%
  na.omit() %>%
  group_by(name) %>%
  summarise(avg = mean(avg)) %>%
  relocate(name) %>%
  filter(name != "")


write.table(pseudobulk, file = "df3_MAGIC_pCb1.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#2
ent <- data.frame(samp2(pseudo, 100, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 700, 23, "Proliferative cell")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 200, 23, "Goblet cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, prol, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

#need to rename genes from ensembl to gene names
ens <- pseudobulk$GENEID

mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
ens_to_symbol_biomart <- getBM(
  filters = 'ensembl_gene_id',
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  values = ens,
  mart = mart)

colnames(ens_to_symbol_biomart) <- c("GENEID", "name")

pseudobulk <- left_join(pseudobulk,ens_to_symbol_biomart, by="GENEID")

#removing any NAs, averaging duplicated genes
#for some reason there is a blank gene named "", I will remove this in the last line of the below code
pseudobulk <- pseudobulk %>%
  na.omit() %>%
  group_by(name) %>%
  summarise(avg = mean(avg)) %>%
  relocate(name) %>%
  filter(name != "")



write.table(pseudobulk, file = "df3_MAGIC_pCb2.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#3
ent <- data.frame(samp2(pseudo, 250, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 550, 23, "Proliferative cell")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 200, 23, "Goblet cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, prol, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"



#need to rename genes from ensembl to gene names
ens <- pseudobulk$GENEID

mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
ens_to_symbol_biomart <- getBM(
  filters = 'ensembl_gene_id',
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  values = ens,
  mart = mart)

colnames(ens_to_symbol_biomart) <- c("GENEID", "name")

pseudobulk <- left_join(pseudobulk,ens_to_symbol_biomart, by="GENEID")

#removing any NAs, averaging duplicated genes
#for some reason there is a blank gene named "", I will remove this in the last line of the below code
pseudobulk <- pseudobulk %>%
  na.omit() %>%
  group_by(name) %>%
  summarise(avg = mean(avg)) %>%
  relocate(name) %>%
  filter(name != "")



write.table(pseudobulk, file = "df3_MAGIC_pCb3.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#4
ent <- data.frame(samp2(pseudo, 750, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 150, 23, "Proliferative cell")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 100, 23, "Goblet cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, prol, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

#need to rename genes from ensembl to gene names
ens <- pseudobulk$GENEID

mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
ens_to_symbol_biomart <- getBM(
  filters = 'ensembl_gene_id',
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  values = ens,
  mart = mart)

colnames(ens_to_symbol_biomart) <- c("GENEID", "name")

pseudobulk <- left_join(pseudobulk,ens_to_symbol_biomart, by="GENEID")

#removing any NAs, averaging duplicated genes
#for some reason there is a blank gene named "", I will remove this in the last line of the below code
pseudobulk <- pseudobulk %>%
  na.omit() %>%
  group_by(name) %>%
  summarise(avg = mean(avg)) %>%
  relocate(name) %>%
  filter(name != "")


write.table(pseudobulk, file = "df3_MAGIC_pCb4.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


#5
ent <- data.frame(samp2(pseudo, 900, 23, "Enterocyte")@assays$MAGIC$data)
prol <- data.frame(samp2(pseudo, 50, 23, "Proliferative cell")@assays$MAGIC$data)
gob <- data.frame(samp2(pseudo, 50, 23, "Goblet cell")@assays$MAGIC$data)

pseudobulk <- cbind(ent, prol, gob) %>%
  mutate_all(function(x) as.numeric((x))) %>%
  rowMeans() %>%
  as.data.frame() %>%
  rownames_to_column(var = "GENEID")
colnames(pseudobulk)[colnames(pseudobulk) == "."] <- "avg"

#need to rename genes from ensembl to gene names
ens <- pseudobulk$GENEID

mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
ens_to_symbol_biomart <- getBM(
  filters = 'ensembl_gene_id',
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  values = ens,
  mart = mart)

colnames(ens_to_symbol_biomart) <- c("GENEID", "name")

pseudobulk <- left_join(pseudobulk,ens_to_symbol_biomart, by="GENEID")

#removing any NAs, averaging duplicated genes
#for some reason there is a blank gene named "", I will remove this in the last line of the below code
pseudobulk <- pseudobulk %>%
  na.omit() %>%
  group_by(name) %>%
  summarise(avg = mean(avg)) %>%
  relocate(name) %>%
  filter(name != "")


write.table(pseudobulk, file = "df3_MAGIC_pCb5.txt", sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


a <- tibble(cellID = colnames(df_samp), clusterID = Idents(df_samp))
b <- data.frame(df_samp@assays$MAGIC$data)
b <- b %>%  mutate_all(function(x) as.numeric((x)))

#b <- expm1(b)
#b[c(1:10), c(1:10)]


#need to rename genes from ensembl to gene names
ens <- rownames(b)
b$GENEID <- rownames(b)

mart <- useDataset('hsapiens_gene_ensembl', useMart('ensembl'))
ens_to_symbol_biomart <- getBM(
  filters = 'ensembl_gene_id',
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  values = ens,
  mart = mart)

colnames(ens_to_symbol_biomart) <- c("GENEID", "name")

b <- left_join(b,ens_to_symbol_biomart, by="GENEID")

#removing any NAs, averaging duplicated genes later below
#for some reason there is a blank gene named "", I will remove this in the last line of the below code
b <- b %>%
  na.omit() %>%
  filter(name != "")%>%
  column_to_rownames(.,"name")%>%
  dplyr::select(.,-GENEID)


##

#b[c(1:10), c(1:10)]
t_a <- data.frame(t(a))
colnames(t_a) <- t_a[1,]

colnames(b) <-gsub(".", "-", colnames(b), fixed = TRUE)


cell_ref <- rbind(t_a, b)


table(colnames(t_a) == colnames(b))
table(colnames(cell_ref) == cell_ref["cellID",])

cell_ref$GENEID <- rownames(cell_ref)
cell_ref <- cell_ref %>%relocate(GENEID)

#Getting rid of the cellID row so the cluster ID is now at the top
cell_ref <- cell_ref[-1,]
##

write.table(cell_ref, 
            file = "Paper3_MAGIC_CellRef_final.txt", 
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
