## cibersort load testis cross compares

if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readxl, ggplot2, ggsci, ggthemes, vroom, readr, tidyverse, DWLS, stringr, Metrics)

fp <- ("~/OneDrive - UW/data for deconv/Deconvolution output (F)/Testis/cibersort_cross_outputs")
files <- list.files(path = fp, pattern = "*.csv", full.names = TRUE)

cb_testis <- list()

for (file in files) {
  data <- read_csv(file)

  basename_file <- tools::file_path_sans_ext(basename(file))
  method_value <- str_extract(basename_file, "(?<=_)[^_]+(?=_testis)")
  basename_minus_method <- str_remove(basename_file, paste0("_", method_value, "\\.csv$"))
  
  data <- data %>%
    { if ("Stroma/Leydig" %in% colnames(.)) rename(., Leydig = `Stroma/Leydig`) else . } %>%
    select(Mixture, Leydig, Sertoli) %>%
    mutate(name = basename_file,
      imp = method_value,
      title = basename_minus_method,
      pb = str_sub(Mixture, -2, -1),
      PBbase = str_extract(Mixture, "^[^_]+"),
      REFbase = str_extract(title, "^[^_]+"),
      cell = ifelse(PBbase == "adult112", Leydig, Sertoli) )%>%
      filter(.,PBbase != 'm_n_d7_SL')
  
  
  cb_testis[[basename_file]] <- data
}

cb_testis <- bind_rows(cb_testis)%>% filter(PBbase != 'm_n_d7_SL')



cb_testis$PBbase<- sub("^(.*_[^_]+)_.*$", "\\1", cb_testis$Mixture)
cb_testis$PBbase <- sub("_[^_]*$", "", cb_testis$PBbase)


cb_testis <- cb_testis %>%
 mutate(PBbase = recode(PBbase,
                "adult112" = "h_a_112", 
                
                "H_i_120" = "h_i_120", 
                "H_I.120" = "h_i_120", 
                "nh120" = "h_i_120", 
                
                "H_i_124" = "h_i_124", 
                "H_I.124" = "h_i_124", 
                "h_n_124" = "h_i_124",
                "hi_124" = "h_i_124",
                
                "m_n_d7_d2" = "m_n_d7_SL",
                
                "m_neo_d2" = "m_n_d2" ),
        REFbase = recode(REFbase,
                         'ha112' = 'h_a_112',
                         'hi120' = 'h_i_120',
                         'hi124' = 'h_i_124',
                         'md2' = 'm_n_d2',
                         'md7' = 'm_d7_diff',
                         'md7' = 'm_n_d7' ))%>%
    filter(PBbase != 'm_n_d7_SL') %>%
  filter(PBbase != 'm_n_d7') %>%
  filter(PBbase != 'm_neo_d7') %>%
  filter(PBbase != REFbase)
                     
saveRDS(cb_testis, "~/OneDrive - UW/data for deconv/Deconvolution output (F)/Testis/cibersort_cross_outputs/cb_testis.rds")
