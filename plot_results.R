#!/usr/bin/env Rscript
# plot_results.R
# Multi-panel summary figure for environmental-filtering workflow outputs.
#
# Usage:
#   Rscript projects/environmental-filtering/plot_results.R \
#     --ef_long    results/EF/EF_long_*.csv \
#     --nst_pair   results/NST/Stochasticity-Ratios_*_PAIRWISE.csv \
#     --assembly   results/QPE/QPE_assembly_processes_*.csv \
#     --output     results/summary_plot.pdf

local({
  conda_prefix <- Sys.getenv("CONDA_PREFIX")
  lib <- if (nchar(conda_prefix) > 0)
    file.path(conda_prefix, "lib", "R", "library")
  else
    .libPaths()[1]
  pkgs <- c("ggplot2", "patchwork")
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) {
    message("Installing into ", lib, ": ", paste(miss, collapse = ", "))
    install.packages(miss, lib = lib, repos = "https://cloud.r-project.org",
                     quiet = TRUE)
  }
})

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--ef_long",
              type = "character", default = NULL,
              help = "Path to EF_long_*.csv (NRI/NTI per sample, long format)"),
  make_option("--nst_pair",
              type = "character", default = NULL,
              help = "Path to Stochasticity-Ratios_*_PAIRWISE.csv"),
  make_option("--assembly",
              type = "character", default = NULL,
              help = "Path to QPE_assembly_processes_*.csv"),
  make_option("--output",
              type = "character", default = "summary_plot.pdf",
              help = "Output path (.pdf or .png) [default: summary_plot.pdf]"),
  make_option("--width",
              type = "double", default = 16,
              help = "Figure width in inches [default: 16]"),
  make_option("--height",
              type = "double", default = 11,
              help = "Figure height in inches [default: 11]"),
  make_option("--group_order",
              type = "character", default = "",
              help = "Comma-separated group names for x-axis ordering (default: alphabetical)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (all(sapply(list(opt$ef_long, opt$nst_pair, opt$assembly), is.null)))
  stop("At least one of --ef_long, --nst_pair, or --assembly must be supplied.")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
resolve_order <- function(groups, opt_order) {
  if (nchar(opt_order) > 0)
    trimws(strsplit(opt_order, ",")[[1]])
  else
    sort(unique(as.character(groups)))
}

# Shared theme: white background, subtle y-grid, rotated x labels
theme_ef <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold"),
      axis.text.x        = element_text(angle = 40, hjust = 1),
      panel.grid.major.y = element_line(colour = "grey93"),
      legend.position    = "right",
      plot.title         = element_text(face = "bold", size = base_size)
    )
}

# Colorblind-friendly palette (up to 12 groups via manual cycling)
group_palette <- c(
  "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
  "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
  "#D37295","#499894"
)

assembly_colours <- c(
  "Homogeneous Selection"   = "#D73027",
  "Variable Selection"      = "#FC8D59",
  "Homogenizing Dispersal"  = "#4575B4",
  "Dispersal Limitation"    = "#91BFDB",
  "Undominated"             = "#CCCCCC"
)

assembly_levels <- names(assembly_colours)

# ---------------------------------------------------------------------------
# Panel A: NRI / NTI
# ---------------------------------------------------------------------------
plot_ef <- NULL

if (!is.null(opt$ef_long) && file.exists(opt$ef_long)) {
  df <- read.csv(opt$ef_long, header = TRUE, stringsAsFactors = FALSE)
  req <- c("value", "measure", "Groups")

  if (!all(req %in% colnames(df))) {
    warning("--ef_long: missing columns (", paste(req, collapse = ", "),
            "). Skipping NRI/NTI panel.")
  } else {
    grp_order  <- resolve_order(df$Groups, opt$group_order)
    df$Groups  <- factor(df$Groups,  levels = grp_order)
    df$measure <- factor(df$measure, levels = c("NRI", "NTI"))
    df <- df[!is.na(df$value) & !is.na(df$Groups), ]

    pal <- setNames(group_palette[seq_along(grp_order)], grp_order)

    plot_ef <- ggplot(df, aes(x = Groups, y = value, fill = Groups)) +
      geom_hline(yintercept = c(-2, 2), linetype = "dashed",
                 colour = "grey40", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "solid",
                 colour = "grey70", linewidth = 0.3) +
      geom_violin(alpha = 0.35, trim = TRUE, colour = NA,
                  data = ~ .x[ave(.x$value, .x$Groups, .x$measure,
                                  FUN = function(v) sum(!is.na(v))) >= 5, ]) +
      geom_boxplot(width = 0.18, outlier.size = 0.6, outlier.alpha = 0.4,
                   colour = "grey25", alpha = 0.8) +
      facet_wrap(~measure, nrow = 1) +
      scale_fill_manual(values = pal, guide = "none") +
      labs(title = "A   Phylogenetic community structure (NRI / NTI)",
           x = NULL, y = "Index value") +
      theme_ef()
  }
}

# ---------------------------------------------------------------------------
# Panel B: NST (pairwise per-sample values)
# ---------------------------------------------------------------------------
plot_nst <- NULL

if (!is.null(opt$nst_pair) && file.exists(opt$nst_pair)) {
  df <- read.csv(opt$nst_pair, header = TRUE, stringsAsFactors = FALSE)

  # Group column is "group" (lowercase) in NST package output
  grp_col <- if ("group" %in% colnames(df)) "group"
              else if ("Groups" %in% colnames(df)) "Groups"
              else NULL

  # NST column: first column matching ^NST (e.g. NST.ij.cao)
  nst_col <- grep("^NST", colnames(df), value = TRUE, ignore.case = FALSE)[1]

  if (is.null(grp_col) || is.na(nst_col)) {
    warning("--nst_pair: could not identify group or NST column. Skipping NST panel.")
  } else {
    df$Groups_ <- as.character(df[[grp_col]])
    df$NST_    <- as.numeric(df[[nst_col]])
    df <- df[!is.na(df$Groups_) & !is.na(df$NST_), ]

    grp_order <- resolve_order(df$Groups_, opt$group_order)
    df$Groups_ <- factor(df$Groups_, levels = grp_order)

    pal <- setNames(group_palette[seq_along(grp_order)], grp_order)

    plot_nst <- ggplot(df, aes(x = Groups_, y = NST_, fill = Groups_)) +
      geom_hline(yintercept = 0.5, linetype = "dashed",
                 colour = "grey40", linewidth = 0.4) +
      geom_violin(alpha = 0.35, trim = TRUE, colour = NA,
                  data = ~ .x[ave(.x$NST_, .x$Groups_,
                                  FUN = function(v) sum(!is.na(v))) >= 5, ]) +
      geom_boxplot(width = 0.18, outlier.size = 0.6, outlier.alpha = 0.4,
                   colour = "grey25", fill = "white", alpha = 0.7) +
      scale_fill_manual(values = pal, guide = "none") +
      scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
      labs(title = "B   Stochasticity (NST)",
           x = NULL, y = "NST") +
      theme_ef()
  }
}

# ---------------------------------------------------------------------------
# Panel C: QPE assembly processes
# ---------------------------------------------------------------------------
plot_qpe <- NULL

if (!is.null(opt$assembly) && file.exists(opt$assembly)) {
  df <- read.csv(opt$assembly, header = TRUE, stringsAsFactors = FALSE)
  req <- c("variable", "percentage", "Groups")

  if (!all(req %in% colnames(df))) {
    warning("--assembly: missing columns (", paste(req, collapse = ", "),
            "). Skipping QPE panel.")
  } else {
    grp_order <- resolve_order(df$Groups, opt$group_order)
    df$Groups   <- factor(df$Groups,   levels = grp_order)
    present_lvls <- assembly_levels[assembly_levels %in% df$variable]
    df$variable <- factor(df$variable, levels = present_lvls)

    plot_qpe <- ggplot(df, aes(x = Groups, y = percentage, fill = variable)) +
      geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
      scale_fill_manual(values = assembly_colours, name = NULL,
                        drop = FALSE, guide = guide_legend(reverse = TRUE)) +
      scale_y_continuous(labels = function(x) paste0(x, "%"),
                         expand = expansion(mult = c(0, 0.02))) +
      labs(title = "C   Community assembly processes (QPE)",
           x = NULL, y = "% of pairwise comparisons") +
      theme_ef() +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.45, "cm"),
            legend.text     = element_text(size = 9))
  }
}

# ---------------------------------------------------------------------------
# Combine panels
# ---------------------------------------------------------------------------
top_panels <- Filter(Negate(is.null), list(plot_ef, plot_nst))

if (length(top_panels) == 0 && is.null(plot_qpe))
  stop("No panels were built. Check that input files exist and have the expected columns.")

combined <- if (length(top_panels) > 0 && !is.null(plot_qpe)) {
  Reduce(`|`, top_panels) / plot_qpe +
    plot_layout(heights = c(1, 1.1))
} else if (length(top_panels) > 0) {
  Reduce(`|`, top_panels)
} else {
  plot_qpe
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
ext <- tolower(tools::file_ext(opt$output))
if (ext == "png") {
  png(opt$output, width = opt$width, height = opt$height,
      units = "in", res = 300)
} else {
  pdf(opt$output, width = opt$width, height = opt$height)
}
print(combined)
dev.off()
message("Saved: ", opt$output)
