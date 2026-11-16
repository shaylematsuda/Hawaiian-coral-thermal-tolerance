# Called by ACE21_16s_analysis.Rmd after constructing the T0 bacterial panels.
# Panel A is exported by ITS2_analysis_update.Rmd; do not rebuild ITS2 data here.
assemble_figure3 <- function(class_plot, genus_plot,
                             its2_export = "Figures/ITS2_Fig3A_plot.rds",
                             output_prefix = "Figures/Fig3") {
  if (!file.exists(its2_export)) {
    stop("Missing ITS2 panel export: ", its2_export,
         ". Run the Figure 3A chunk in ITS2_analysis_update.Rmd first.")
  }
  panel_a <- readRDS(its2_export)
  stopifnot(inherits(panel_a, "ggplot"))
  # Draw each displayed category once per host. Stacking many pooled OTU
  # rectangles creates fine PDF rendering seams inside the grey category.
  solid_categories <- function(plot) {
    original_totals <- tapply(plot$data$mean, plot$data$Species, sum)
    plot$data <- plot$data |>
      dplyr::group_by(Species, display_taxon) |>
      dplyr::summarise(mean = sum(mean), .groups = "drop")
    stopifnot(isTRUE(all.equal(
      original_totals, tapply(plot$data$mean, plot$data$Species, sum))))
    plot
  }
  class_plot <- solid_categories(class_plot)
  genus_plot <- solid_categories(genus_plot)
  # Panel A is already assembled in the ITS2 code (DIVs, profiles, ALG labels).
  # Preserve its component layout and profile-only legend here.
  bacterial_theme <- ggplot2::theme(
    text = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(size = 9),
    axis.title.y = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 7),
    legend.title = ggplot2::element_text(size = 10),
    legend.key.height = grid::unit(3.5, "mm"),
    legend.key.width = grid::unit(3.5, "mm"),
    plot.margin = ggplot2::margin(12, 8, 8, 12))
  # Retain the original genus-panel palette and use readable legend labels.
  genus_colors <- c(
    "<2%" = "#AAAAAB", "Alphaproteobacteria_unclassified" = "#2B8CBC",
    "Bacteria_unclassified" = "#9ABCAD", "Candidatus_Amoebophilus" = "#D6BAD9",
    "Cytophagales_unclassified" = "#F6BD91", "Endozoicomonas" = "#9FD7F5",
    "Fulvivirga" = "#AA76B3", "Gammaproteobacteria_unclassified" = "#70BADD",
    "Mastigocoleus_BC008" = "#AA7552", "MBIC10086" = "#C97D67",
    "P3OB-42_ge" = "#8D6373", "Rhodobacteraceae_unclassified" = "#DBF1F9",
    "uncultured" = "#476269")
  genus_taxa <- sort(unique(as.character(genus_plot$data$display_taxon)))
  extra_taxa <- setdiff(genus_taxa, names(genus_colors))
  if (length(extra_taxa)) {
    genus_colors <- c(genus_colors, stats::setNames(
      grDevices::hcl.colors(length(extra_taxa), "Dark 3"), extra_taxa))
  }
  genus_plot <- genus_plot + ggplot2::scale_fill_manual(
    values = genus_colors,
    labels = function(x) {
      x <- sub("_unclassified$", " (unclassified)", x)
      gsub("_", " ", x)
    })
  bacterial_y <- ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, .25),
    expand = ggplot2::expansion(mult = c(0, .01)))
  # Formatting only: use wider bars with visible gaps; preserve x-axis extents.
  narrow_bars <- function(plot) {
    plot$layers <- lapply(plot$layers, function(layer) {
      new_layer <- ggplot2::ggproto(NULL, layer)
      if (inherits(layer$geom, "GeomCol") || inherits(layer$geom, "GeomBar")) {
        new_layer$geom_params <- modifyList(layer$geom_params, list(width = .7))
        # Match the ITS2 profile outlines; categories are already pooled above.
        new_layer$aes_params <- modifyList(layer$aes_params,
                                           list(colour = "black", linewidth = .08))
      }
      new_layer
    })
    plot
  }
  class_plot <- narrow_bars(class_plot)
  genus_plot <- narrow_bars(genus_plot)
  # Reference layout: one B row with a right-hand legend beside each chart.
  compact_theme <- ggplot2::theme(
    text = ggplot2::element_text(size = 9),
    axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 8),
    axis.title.y = ggplot2::element_text(size = 9),
    legend.position = "right", legend.text = ggplot2::element_text(size = 5.5),
    legend.title = ggplot2::element_text(size = 10),
    legend.key.height = grid::unit(1.7, "mm"),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.spacing.y = grid::unit(0, "mm"),
    legend.spacing.y = grid::unit(0, "mm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    plot.margin = ggplot2::margin(6, 2, 6, 2))
  class_panel <- suppressMessages(class_plot + bacterial_y + compact_theme) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 7),
                   legend.key.height = grid::unit(5, "mm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
  genus_panel <- suppressMessages(genus_plot + bacterial_y + compact_theme) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 6.5),
                   legend.key.height = grid::unit(5, "mm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
  right_legend <- function(plot) {
    grob <- ggplot2::ggplotGrob(plot)
    index <- which(grob$layout$name == "guide-box-right")
    stopifnot(length(index) == 1)
    grob$grobs[[index]]
  }
  class_legend <- right_legend(class_panel)
  genus_legend <- right_legend(genus_panel)
  charts <- cowplot::align_plots(
    class_panel + ggplot2::theme(legend.position = "none"),
    genus_panel + ggplot2::theme(legend.position = "none"),
    align = "hv", axis = "tblr")
  bottom <- cowplot::plot_grid(
    charts[[1]], class_legend, charts[[2]], genus_legend,
    nrow = 1, rel_widths = c(1, .88, 1, 1.27))
  figure <- cowplot::plot_grid(
    panel_a, bottom, ncol = 1, rel_heights = c(1.45, 1),
    labels = c("A", "B"), label_size = 15,
    hjust = 0, label_x = 0, label_y = 1)
  ggplot2::ggsave(paste0(output_prefix, ".pdf"), figure,
                  width = 8, height = 10, bg = "white")
  ggplot2::ggsave(paste0(output_prefix, ".png"), figure,
                  width = 8, height = 10, dpi = 300, bg = "white")
  saveRDS(figure, paste0(output_prefix, "_plot.rds"))
  figure
}
