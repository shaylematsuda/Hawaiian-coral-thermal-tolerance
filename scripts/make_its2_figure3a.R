# Shared plotting function used by the Figure 3A chunk in the ITS2 notebook.
make_its2_figure3a <- function(profile_ps, div_ps) {
  species_levels <- c("Montipora_capitata", "Porites_compressa",
                      "Pavona_varians", "Pocillopora_acuta")
  species_labels <- c("Montipora capitata", "Porites compressa", "Pavona varians", "Pocillopora acuta")
  alg_levels <- c("C31", "C31D4", "C15cj", "C15n", "C1", "C27", "C1d", "C42_2")
  profiles <- phyloseq::subset_samples(profile_ps, Time.Point == "T0")
  divs <- phyloseq::subset_samples(div_ps, Time.Point == "T0")
  # Missing samples must be resolved explicitly rather than silently dropping bars.
  stopifnot(setequal(phyloseq::sample_names(profiles), phyloseq::sample_names(divs)))
  meta <- as(phyloseq::sample_data(profiles), "data.frame")
  meta$Sample <- rownames(meta)
  meta$Species <- factor(meta$Species, levels = species_levels)
  meta$ALG <- factor(as.character(meta$Clade), levels = alg_levels)
  stopifnot(!anyNA(meta$ALG), !anyNA(meta$Species), !anyDuplicated(meta$Sample))
  meta <- meta[order(meta$Species, meta$ALG, meta$Sample), ]
  meta$x <- ave(seq_len(nrow(meta)), meta$Species, FUN = seq_along)
  meta$ALG_label <- sub("C42_2", "C42.2", as.character(meta$ALG))

  melt_rank <- function(ps, rank) {
    ps <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
    long <- phyloseq::psmelt(ps)
    long$Sample <- as.character(long$Sample)
    long$x <- meta$x[match(long$Sample, meta$Sample)]
    long$Species <- factor(long$Species, levels = species_levels)
    totals <- tapply(long$Abundance, long$Sample, sum)
    stopifnot(all(abs(totals - 1) < 1e-10), !anyNA(long$x))
    long
  }
  profile_data <- melt_rank(profiles, "its2_type_profile")
  div_data <- melt_rank(divs, "DIV")
  profile_names <- sort(unique(as.character(profile_data$its2_type_profile)))
  c_names <- profile_names[!grepl("^D", profile_names)]
  d_names <- profile_names[grepl("^D", profile_names)]
  profile_colors <- c(
    stats::setNames(grDevices::colorRampPalette(c("lightblue", "blue", "darkblue"))(length(c_names)), c_names),
    stats::setNames(grDevices::colorRampPalette(c("lightyellow", "orange", "darkorange"))(length(d_names)), d_names))
  div_names <- sort(unique(as.character(div_data$DIV)))
  palette <- scales::hue_pal()(length(div_names))
  div_colors <- stats::setNames(palette, div_names)
  facet <- ggplot2::facet_wrap(~Species, nrow = 1, scales = "free_x",
    labeller = ggplot2::as_labeller(stats::setNames(species_labels, species_levels)))
  common <- ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold.italic", size = 10),
      panel.spacing.x = grid::unit(2, "mm"),
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 0, 3))
  xscale <- ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = .5))
  div_plot <- ggplot2::ggplot(div_data, ggplot2::aes(x, Abundance, fill = DIV)) +
    ggplot2::geom_col(width = .95, color = "black", linewidth = .08) + facet + xscale +
    ggplot2::scale_fill_manual(values = div_colors, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
      breaks = seq(0, 1, .25), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = NULL, y = "Relative Abundance") + common +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 3, 6, 3))
  profile_plot <- ggplot2::ggplot(profile_data,
      ggplot2::aes(x, Abundance, fill = its2_type_profile)) +
    ggplot2::geom_col(width = .95, color = "black", linewidth = .08) + facet + xscale +
    ggplot2::scale_fill_manual(values = profile_colors, name = "Algal Symbiont Lineages") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
      breaks = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = NULL, y = NULL) + common +
    ggplot2::theme(strip.text = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 6, 3),
      legend.text = ggplot2::element_text(size = 5.2),
      legend.title = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      legend.key.spacing.y = grid::unit(0, "mm"),
      legend.spacing.y = grid::unit(0, "mm"),
      legend.key.height = grid::unit(2.8, "mm"),
      legend.key.width = grid::unit(3, "mm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3, title.position = "top"))
  legend_grob <- ggplot2::ggplotGrob(profile_plot)
  legend <- legend_grob$grobs[[which(legend_grob$layout$name == "guide-box-bottom")]]
  profile_plot <- profile_plot + ggplot2::theme(legend.position = "none")
  groups <- meta |>
    dplyr::group_by(Species, ALG_label) |>
    dplyr::summarise(left = min(x) - .4, right = max(x) + .4,
                     middle = mean(range(x)), .groups = "drop")
  annotation <- ggplot2::ggplot(meta, ggplot2::aes(x = x)) + facet + xscale +
    ggplot2::geom_segment(data = groups,
      ggplot2::aes(x = left, xend = right, y = .65, yend = .65),
      inherit.aes = FALSE, linewidth = .5) +
    ggplot2::geom_text(data = groups,
      ggplot2::aes(x = middle, y = .22, label = ALG_label),
      inherit.aes = FALSE, size = 2.7) +
    ggplot2::scale_y_continuous(limits = c(-.05, 1.05), expand = c(0, 0)) +
    common + ggplot2::theme(
      strip.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank())
  rows <- cowplot::align_plots(div_plot, profile_plot, annotation, align = "v", axis = "lr")
  left <- cowplot::plot_grid(plotlist = rows, ncol = 1, rel_heights = c(4, 1, .35))
  panel <- cowplot::plot_grid(left, legend, ncol = 1, rel_heights = c(1, .48))
  list(plot = panel, div_plot = div_plot, profile_plot = profile_plot,
       sample_order = meta, alg_groups = groups,
       profile_data = profile_data, div_data = div_data)
}
