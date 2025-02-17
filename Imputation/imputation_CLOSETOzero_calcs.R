### Get the # of zeros in each matrix with imputation

library(data.table)

#function
getalmostzero <- function(i){
almost_zeros <- sum(i > 0 & i < 100)
total_elements <- length(i)
per_zero <- (almost_zeros / total_elements) * 100

return(per_zero)
}


#no imputation
h_120 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/objects/h_infant_120_data.rds")
h_120_zero <- getalmostzero(h_120)

h_124 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/objects/h_infant_124_data.rds")
h_124_zero <- getalmostzero(h_124)

h_112 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/objects/h_adult_112_data.rds")
h_112_zero <- getalmostzero(h_112)

mn2 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/objects/m_neonate_d2_data.rds")
mn2_zero <- getalmostzero(mn2)

mn7 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/objects/m_neonate_d7_data.rds")
mn7_zero <- getalmostzero(mn7)


##ALRA
h_120 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/alra/h_infant_120_t_alra.rds")[[3]]
h_120_zero_alra <- getalmostzero(h_120)

h_124 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/alra/h_infant_124_t_alra.rds")[[3]]
h_124_zero_alra <- getalmostzero(h_124)

h_112 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/alra/h_adult_112_t_alra.rds")[[3]]
h_112_zero_alra <- getalmostzero(h_112)

mn2 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/alra/m_neonate_d2_t_alra.rds")[[3]]
mn2_zero_alra <- getalmostzero(mn2)

mn7 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/alra/m_neonate_d7_t_alra.rds")[[3]]
mn7_zero_alra <- getalmostzero(mn7)


##SAVER
h_120 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/saver/hi120_saver_cpm_SAVERoutput.rds")[[1]]
h_120_zero_saver <- getalmostzero(h_120)

h_124 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/saver/hi124_saver_cpm_SAVERoutput.rds")[[1]]
h_124_zero_saver <- getalmostzero(h_124)

h_112 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/saver/ha112_saver_cpm_SAVERoutput.rds")[[1]]
h_112_zero_saver <- getalmostzero(h_112)

mn2 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/saver/md2_saver_cpm_SAVERoutput.rds")[[1]]
mn2_zero_saver <- getalmostzero(mn2)

mn7 <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/saver/md7_saver_cpm_SAVERoutput.rds")[[1]]
mn7_zero_saver <- getalmostzero(mn7)

##MAGIC
h_120 <- as.matrix(fread("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/magic/hi120_MAGIC_matrix.csv"))
h_120_zero_magic <- getalmostzero(h_120)
rm(h_120)


h_124 <- as.matrix(fread("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/magic/hi124_MAGIC_matrix.csv"))
h_124_zero_magic <- getalmostzero(h_124)
rm(h_124)

h_112 <- fread("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/magic/ha112_MAGIC_matrix.csv")
h_112_zero_magic <- getalmostzero(h_112)

rm(h_112)

mn2 <- fread("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/magic/mn2_MAGIC_matrix.csv")
mn2_zero_magic <- getalmostzero(mn2)

rm(mn2)

mn7 <- fread("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/magic/mn7_MAGIC_matrix.csv")
mn7_zero_magic <- getalmostzero(mn7)


