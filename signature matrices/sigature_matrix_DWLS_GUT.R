## DWLS Signature Matrix

library(DWLS)
library(dplyr)

# List of input files (you can adjust the path or list method)
file_list <- list.files(path = "~/OneDrive - UW/data for deconv/Paper revision materials/DATA/cpm/GUT/Gut_cellRefs", pattern = "*MAGIC_CellRef_final.txt", full.names = TRUE, recursive = FALSE)

for (file in file_list) {
  
  df <- read.delim(file_list[1], header=FALSE)
  
  rownames(df) <- df$V1
  df <- df[,-1]
  
  ids <- df[1,]
  df <- df[-1,]
  
  df <- df %>% mutate_all(function(x) as.numeric(x))
  
  ids <- gsub(" ", "", ids, fixed = TRUE)
  ids <- gsub("/", "", ids, fixed = TRUE)
  ids <- gsub("-", "", ids, fixed = TRUE)
  
  output_dir <- file.path(dirname(file), tools::file_path_sans_ext(basename(file)))
  dir.create(output_dir, showWarnings = FALSE)
  
  Signature <- buildSignatureMatrixMAST(df,  ids,  path = output_dir,  diff.cutoff = 0.5, pval.cutoff = 0.01, f = 300)
  
  
  
}




