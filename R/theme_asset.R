# ==============================================================
# R/theme_asset.R — shared plotting theme and palettes
#
# One theme, applied everywhere, so figures are internally consistent
# and journal-ready without per-figure fiddling.
# ==============================================================

#' Publication theme
#'
#' Sized for Arthritis & Rheumatology / JCI Insight single-column
#' figures: 7pt base, thin lines, no chartjunk.
theme_asset <- function(base_size = 7, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2,
                                               colour = "grey92"),
      panel.border     = ggplot2::element_rect(linewidth = 0.4,
                                               colour = "grey20"),
      axis.ticks       = ggplot2::element_line(linewidth = 0.3),
      axis.text        = ggplot2::element_text(colour = "grey20"),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold",
                                               size = base_size),
      legend.key.size  = grid::unit(0.35, "cm"),
      legend.background = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold",
                                               size = base_size + 1),
      plot.subtitle    = ggplot2::element_text(colour = "grey35"),
      plot.caption     = ggplot2::element_text(colour = "grey45",
                                               size = base_size - 1,
                                               hjust = 0)
    )
}

# ---- palettes ------------------------------------------------
# Fixed colours so the same category is the same colour in every figure.

PAL_ARM <- c(PLACEBO = "#3B6EA5", ABATACEPT = "#C0504D", HC = "#7F7F7F")

PAL_GROUP <- c(SSC = "#C0504D", HC = "#7F7F7F")

PAL_TIME <- c(M00 = "#D9D9D9", M03 = "#9ECAE1", M06 = "#2B7BBA")

PAL_IMPROVER <- c(improver = "#2A9D8F", non_improver = "#E76F51")

# Intrinsic subsets — matching the convention in the Whitfield-lab
# figures so panels read the same way across papers.
PAL_SUBSET <- c(inflammatory      = "#6A3D9A",
                fibroproliferative = "#E31A1C",
                normal_like        = "#33A02C",
                limited            = "#FF7F00")

PAL_LINEAGE <- c(Stromal      = "#E31A1C",
                 Myeloid      = "#FF7F00",
                 TNK          = "#1F78B4",
                 BPlasma      = "#6A3D9A",
                 Vascular     = "#33A02C",
                 Keratinocyte = "#B15928",
                 Other        = "#BDBDBD")

#' Colourblind-safe qualitative palette for arbitrary categories
pal_asset <- function(n) {
  base <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
            "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD")
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::colorRampPalette(base)(n)
}

# ---- figure sizes --------------------------------------------
# Journal column widths in inches. Use these rather than guessing.
FIG_W_SINGLE <- 3.42   # 87 mm
FIG_W_DOUBLE <- 7.09   # 180 mm

#' Save a figure at a journal column width
#'
#' @param p ggplot
#' @param name file stem, e.g. "fig2c_fibroblast_composition_by_improver"
#' @param width one of "single", "double", or a number in inches
save_figure <- function(p, name, width = "single", height = 3,
                        dir = FIGURES, device = "pdf") {
  w <- switch(as.character(width),
              single = FIG_W_SINGLE,
              double = FIG_W_DOUBLE,
              as.numeric(width))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, paste0(name, ".", device))
  ggplot2::ggsave(f, p, width = w, height = height, device = device,
                  useDingbats = FALSE)
  message("  figure -> ", f)
  invisible(f)
}
