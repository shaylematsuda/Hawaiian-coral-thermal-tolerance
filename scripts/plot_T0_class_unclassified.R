# Regenerate only the T0 class composition using the current notebook functions.
# Run from the repository root: Rscript scripts/plot_T0_class_unclassified.R
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
})
load("Data/Bac.seq_phyloseq.RData")
excluded <- c("S1753b", "S208", "S3476", "SN35b", "S3149", "S1240")
ps <- subset_samples(Bac.seq, Type == "sample" & Time.Point == "T0" &
                       !sample_name.1 %in% excluded)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- transform_sample_counts(ps, function(x) x / sum(x))

# Use exactly the plotting and summarizing functions maintained in the notebook.
notebook <- readLines("ACE21_16s_analysis.Rmd")
start <- which(notebook == 'h16_species_levels <- c(')
end <- which(notebook == 'h16_ambient_rank_summaries <- lapply(')
stopifnot(length(start) == 1, length(end) == 1, end > start)
eval(parse(text = notebook[start:(end - 1)]))
summary <- h16_summarize_rank(ps, "Class")
plot <- h16_plot_rank(summary, "Class") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01)))

# Every bar must sum to 100%; confirm the previously calculated PVAR category.
totals <- tapply(summary$mean, summary$Species, sum)
stopifnot(all(abs(totals - 1) < 1e-10))
pvar <- subset(summary, Species == "Pavona_varians" & Class == "Bacteria_unclassified")
stopifnot(nrow(pvar) == 1, abs(pvar$mean - 0.0713333333333333) < 1e-10, pvar$N == 12)

prefix <- "Figures/Figure3_T0_class_with_unclassified"
ggsave(paste0(prefix, ".png"), plot, width = 8, height = 4.5, dpi = 300, bg = "white")
ggsave(paste0(prefix, ".pdf"), plot, width = 8, height = 4.5, bg = "white")
write.csv(summary, paste0(prefix, "_data.csv"), row.names = FALSE)
print(subset(summary, Class == "Bacteria_unclassified", c("Species", "N", "mean")))
