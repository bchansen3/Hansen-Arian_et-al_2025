
#load packages

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, gridExtra, lme4, scales, tidyverse, ggfortify,cluster, ggthemes, ggpubr, ggprism, extrafont, readxl, coin, ggsci)

theme_set(theme_few() + theme(text = element_text(size=28, family="Myriad Pro")))

### qPCR gene plots

qpcr <- read.csv("~/OneDrive - UW/research/qpcr/qpcr/fold_changes_4-15_bioreps1-2-3.csv")

# Create the plot with individual points and a line through the means
a <- ggplot(plot_testos %>% filter(gene_name == genes[i]), aes(x = time, y = logcpm)) +
  geom_point(aes(color = rep), size = 2.75, alpha = 0.65, shape =15) +  # Plot individual points with transparency
  geom_line(data = plot_summary, aes(x = time, y = mean_logcpm), color = "black", size = 2, linewidth = 1.2, alpha = 0.9, linetype = "dashed") +  # Line through means
  labs(x = "Time", y = "Testosterone (ng/mL)") +
  scale_x_continuous(breaks = c(0, 2, 6, 12, 48)) +
  theme_classic(base_size = 20, base_family = "Myriad Pro") +
  scale_color_d3() +
  theme(
    legend.position = "none", # Hide the legend
    panel.grid = element_blank() # Remove grid lines for cleaner look
  )

# Create title with gene name in italics
plot_title <- paste0(genes[i])  # Concatenate gene name with title string

# Apply the title to the plot
a <- a + ggtitle(plot_title)  # Add title using the constructed string


ggsave(filename = paste0("~/OneDrive - UW/data for deconv/Deconvolution output (F)/Testis/plots/", genes[i], "_induction.png"), plot = a,dpi = 600, width = 5, height = 4, units = 'in')
}

