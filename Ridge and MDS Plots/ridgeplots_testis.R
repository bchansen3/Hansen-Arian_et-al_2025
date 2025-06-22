### Make RIDGE plots for each gene
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readxl, ggplot2, ggsci, ggthemes, vroom, readr, tidyverse, DWLS, stringr, janitor, ggridges, purrr, Seurat)

# need to load in each reference


fp <- ("~/OneDrive - UW/data for deconv/Paper revision materials/DATA/cpm/Testis/cellrefs")
files <- list.files(path = fp, pattern = "*.txt", full.names = TRUE)
gene_df <- list()

for (file in files) {
  
  basename_file <- basename(file)
  
  ref <- str_remove(basename_file, paste0("_", method_value, "\\.txt$"))
  
  df <- read.delim(file, header = FALSE) %>%
    t() %>%
    row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    select(Sox9, Inhba, Dazl, Ddx4, clusterID) %>%
    mutate(
      method = method_value,
      ref = ref,
      type = clusterID,  
      across(c(Sox9, Inhba, Dazl, Ddx4), as.numeric),  
      across(c(Sox9, Inhba, Dazl, Ddx4), ~ log2(. + 1)))
  
  gene_df[[length(gene_df) + 1]] <- df
}

gene_plots <- bind_rows(gene_df)

gene_plots <- gene_plots%>% mutate(method = sub(".*(ALRA|SAVER|SVR|MAGIC|noIMP).*", "\\1", ref),
                                   method = case_when(method == "noIMP" ~ "No Imputation", TRUE ~ method),
                                   method = case_when(method == "SVR" ~ "SAVER", TRUE ~ method))


gene_plots$method <- factor(gene_plots$method, levels = c("No Imputation", "ALRA", "SAVER", "MAGIC"))



table(gene_plots$method)

input <- gene_plots %>% mutate(type = case_when(
                               type %in% c("d_spg", "Germ", "Pro", 'Pro-SG/SSC',"SSCs") ~ "Germ", 
                               type %in% c("Leydig", "Stroma", "Stroma/Leydig") ~ "Leydig",
                               TRUE ~ type),
                               species = case_when(
                                 grepl("h", ref, ignore.case = TRUE) ~ "Human",
                                 grepl("m", ref, ignore.case = TRUE) ~ "Mouse",
                                 TRUE ~ NA)) %>% 
                                      filter(type %in% c("Germ", "Leydig", "Sertoli"))


genelist <- c('Sox9', 'Inhba', 'Dazl', 'Ddx4')

for (gene in genelist) {
  
  ridge <- ggplot(data=input, aes(x = get(gene), y = type, fill = species)) +
    geom_density_ridges(alpha = 0.65, scale = 1.1) +
    scale_y_discrete(expand = expansion(add = c(0.05, 1.2))) +
    ylab(gene) +  
    xlab("logCPM") +
    coord_cartesian(clip = "off") +
    theme_bw(28, base_family = "Myriad Pro") +
    scale_fill_d3()+
    theme(legend.position = "bottom") +
    theme(axis.title.y=element_text(face="italic"), base_family = "Myriad Pro")+
    theme(legend.title=element_blank())+
    facet_wrap(~method, nrow = 1)
  
  output_file <- paste0("~/OneDrive - UW/data for deconv/Paper revision materials/updated_figures/testis_ridges/testis_", gene, "_ridge.png")
  ggsave(ridge, dpi = 600, filename = output_file, width = 14, height = 7)
  
}

