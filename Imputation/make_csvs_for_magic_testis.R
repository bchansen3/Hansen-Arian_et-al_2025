# use this to make .csv files for MAGIC

h_infant_120_t <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/h_infant_120_data.rds")%>% t()
write.csv(h_infant_120_t, "~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/h_infant_120_data_transposed.csv", row.names = TRUE)


h_infant_124_t <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/h_infant_124_data.rds")%>% t()
write.csv(h_infant_124_t, "~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/h_infant_124_t.csv", row.names = TRUE)

h_adult_112_t <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/h_adult_112_data.rds")%>% t()
write.csv(h_adult_112_t, "~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/h_adult_112_t.csv", row.names = TRUE)

m_neonate_d2_t <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/m_neonate_d2_data.rds")%>% t()
write.csv(m_neonate_d2_t, "~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/m_neonate_d2_t.csv", row.names = TRUE)

m_neonate_d7_t <- readRDS("~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/m_neonate_d7_data.rds")%>% t()
write.csv(m_neonate_d7_t, "~/OneDrive - UW/data for deconv/Data matrices (C)/Testis/matricesforsaver/m_neonate_d7_t.csv", row.names = TRUE)
