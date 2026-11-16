# Rebuild the exported ITS2 panel and Figure 3 without running unrelated analyses.
# Run from the repository root: Rscript scripts/render_figure3.R
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
})
# Execute the existing ITS2 Figure 3A chunk, including its panel export.
coral_profile2 <- readRDS("coral_profile2_cleaned.RDS")
coral_DIV2_rare <- readRDS("coral_DIV2_rare.RDS")
its2_notebook <- readLines("ITS2_analysis_update.Rmd")
start <- which(grepl("^```\\{r figure-3a-divs-and-profiles", its2_notebook)) + 1
end <- which(seq_along(its2_notebook) > start & its2_notebook == "```")[1] - 1
stopifnot(length(start) == 1, length(end) == 1)
eval(parse(text = its2_notebook[start:end]))

# Reuse the filtered T0 data and notebook functions used for the class preview.
source("scripts/plot_T0_class_unclassified.R")
class_plot <- h16_plot_rank(h16_summarize_rank(ps, "Class"), "Class")
genus_plot <- h16_plot_rank(h16_figure3_genus_data(ps), "Genus")
write.csv(class_plot$data, "Figures/Fig3_class_data.csv", row.names = FALSE)
write.csv(genus_plot$data, "Figures/Fig3_genus_data.csv", row.names = FALSE)
stopifnot(all(abs(tapply(genus_plot$data$mean, genus_plot$data$Species, sum) - 1) < 1e-10))
source("scripts/assemble_figure3.R")
Fig3 <- assemble_figure3(class_plot, genus_plot)
