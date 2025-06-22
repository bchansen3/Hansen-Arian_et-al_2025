
library(DWLS)
library(dplyr)

# List of input files (you can adjust the path or list method)
file_list <- list.files(path = "~/OneDrive - UW/data for deconv/Signature matrices impute+no impute (D)/Testis/DWLS_testis/alldata_revision_04132025/", pattern = "*MAGIC.txt", full.names = TRUE)

for (file in file_list) {
  
  df <- read.delim(file, header=FALSE)
  
  rownames(df) <- df$V1
  df <- df[,-1]
  
  ids <- df[1,]
  df <- df[-1,]
  
  df <- df %>% mutate_all(function(x) as.numeric(x))
  
  #df <- log(df + 1, base = 2)
  
  ids <- gsub(" ", "", ids, fixed = TRUE)
  ids <- gsub("/", "", ids, fixed = TRUE)
  ids <- gsub("-", "", ids, fixed = TRUE)
  
  output_dir <- file.path(dirname(file), tools::file_path_sans_ext(basename(file)))
  dir.create(output_dir, showWarnings = FALSE)
  
  Signature <- buildSignatureMatrixMAST(df,  ids,  path = output_dir,  diff.cutoff = 0.5, pval.cutoff = 0.01, f = 300)
  
  
  
}






