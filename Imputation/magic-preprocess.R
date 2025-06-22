#MAGIC change DFs to TSV files

library(Seurat)

df1_ds <- readRDS( "~/OneDrive - UW/data for deconv/CPM SEURAT OBJECTS/DS_paper1_cpm.rds")
df2_ds <- readRDS( "~/OneDrive - UW/data for deconv/CPM SEURAT OBJECTS/DS_paper2_cpm.rds")
df3_ds <- readRDS( "~/OneDrive - UW/data for deconv/CPM SEURAT OBJECTS/DS_paper3_cpm.rds")

df1_ds <- NormalizeData(df1_ds, normalization.method="RC", scale.factor=1e6)
df2_ds <- NormalizeData(df2_ds, normalization.method="RC", scale.factor=1e6)
df3_ds <- NormalizeData(df3_ds, normalization.method="RC", scale.factor=1e6)


df1_matrix <- as.matrix(GetAssayData(object = df1_ds, layer = "data"))
df2_matrix <- as.matrix(GetAssayData(object = df2_ds, layer = "data"))
df3_matrix <- as.matrix(GetAssayData(object = df3_ds, layer = "data"))


write.csv(t(df1_matrix), "~/OneDrive - UW/data for deconv/Data matrices (C)/Gut/dataMAGIC/Tdf1_dsmatrix.csv", row.names = TRUE)
write.csv(t(df2_matrix), "~/OneDrive - UW/data for deconv/Data matrices (C)/Gut/dataMAGIC/Tdf2_dsmatrix.csv", row.names = TRUE)
write.csv(t(df3_matrix), "~/OneDrive - UW/data for deconv/Data matrices (C)/Gut/dataMAGIC/Tdf3_dsmatrix.csv", row.names = TRUE)