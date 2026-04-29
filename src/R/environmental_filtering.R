# environmental_filtering.R
# NRI/NTI phylogenetic alpha diversity — environmental filtering analysis.
# Outputs three CSVs: wide NRI/NTI, long format with metadata, pairwise test results.

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(ape)
  library(picante)
})

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option_list <- list(
  # Input / output
  make_option("--feature_table",
              type="character", default=NULL,
              help="Path to feature table (BIOM, TSV, or GTDB TSV) [required]"),
  make_option("--biom_file",
              type="character", default=NULL,
              help="Alias for --feature_table (deprecated)"),
  make_option("--meta_table",
              type="character", default=NULL,
              help="Path to sample metadata CSV [required]"),
  make_option("--tree_file",
              type="character", default=NULL,
              help="Path to phylogenetic tree in Newick format [required]"),
  make_option("--taxonomy_table",
              type="character", default="",
              help="Path to two-column TSV taxonomy table (required for tsv/gtdb formats)"),
  make_option("--input_format",
              type="character", default="biom",
              help="Input format: biom | tsv | gtdb [default: biom]"),
  make_option("--output_dir",
              type="character", default=".",
              help="Directory for output CSVs [default: .]"),
  make_option("--label",
              type="character", default="analysis",
              help="Label appended to output filenames [default: analysis]"),
  make_option("--scripts_dir",
              type="character", default="/opt/ecology-scripts",
              help="Path to ecology-scripts root inside the container"),

  # Filtering
  make_option("--min_library_size",
              type="integer", default=5000L,
              help="Minimum read depth per sample [default: 5000]"),
  make_option("--exclude_column",
              type="character", default="",
              help="Metadata column used to exclude samples"),
  make_option("--exclude_values",
              type="character", default="",
              help="Comma-separated values in exclude_column to remove"),

  # Grouping
  make_option("--group",
              type="character", default="",
              help="One metadata column, or comma-separated columns to paste (space-joined) as the Groups label"),
  make_option("--type",
              type="character", default="",
              help="Metadata column for sample Type (included in output CSVs)"),
  make_option("--type2",
              type="character", default="",
              help="Metadata column for secondary Type2 annotation (included in output CSVs)"),

  # NRI / NTI
  make_option("--runs",
              type="integer", default=999L,
              help="Number of null model randomizations [default: 999]"),
  make_option("--iterations",
              type="integer", default=1000L,
              help="Iterations for trialswap null model [default: 1000]"),
  make_option("--top_n_features",
              type="integer", default=2000L,
              help="Number of most abundant features to use [default: 2000]"),
  make_option("--null_model",
              type="character", default="trialswap",
              help="Null model: taxa.labels|richness|frequency|sample.pool|phylogeny.pool|independentswap|trialswap [default: trialswap]"),
  make_option("--abundance_weighted",
              action="store_true", default=FALSE,
              help="Use abundance-weighted metrics [default: TRUE via nextflow.config]"),

  # Statistical testing
  make_option("--test_method",
              type="character", default="anova",
              help="Pairwise test: anova | kruskal | none [default: anova]"),
  make_option("--p_adjust_method",
              type="character", default="BH",
              help="P-value adjustment method: BH | bonferroni | holm | none [default: BH]"),
  make_option("--threads",
              type="integer", default=1L,
              help="Number of parallel threads (ses.mpd and ses.mntd run concurrently when >= 2) [default: 1]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# biom_file alias fallback
if (is.null(opt$feature_table) && !is.null(opt$biom_file))
  opt$feature_table <- opt$biom_file

# Validation
# ---------------------------------------------------------------------------
if (is.null(opt$feature_table))
  stop("--feature_table (or --biom_file) is required.")
if (is.null(opt$meta_table))
  stop("--meta_table is required.")
if (is.null(opt$tree_file) || opt$tree_file == "")
  stop("--tree_file is required for NRI/NTI calculation.")
if (!opt$input_format %in% c("biom", "tsv", "gtdb"))
  stop("--input_format must be one of: biom, tsv, gtdb.")
if (!opt$null_model %in% c("taxa.labels","richness","frequency","sample.pool",
                             "phylogeny.pool","independentswap","trialswap"))
  stop("--null_model must be one of: taxa.labels, richness, frequency, sample.pool, ",
       "phylogeny.pool, independentswap, trialswap.")
if (!opt$test_method %in% c("anova", "kruskal", "none"))
  stop("--test_method must be one of: anova, kruskal, none.")
if (!opt$p_adjust_method %in% c("BH", "bonferroni", "holm", "none"))
  stop("--p_adjust_method must be one of: BH, bonferroni, holm, none.")

if (!dir.exists(opt$output_dir))
  dir.create(opt$output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load helper
# ---------------------------------------------------------------------------
source(file.path(opt$scripts_dir, "src", "R", "load_feature_table.R"))

# ---------------------------------------------------------------------------
# Feature table + taxonomy
# ---------------------------------------------------------------------------
message("Loading feature table...")
tax_arg <- if (nchar(opt$taxonomy_table) > 0) opt$taxonomy_table else NULL
ft_obj  <- load_feature_table(opt$feature_table, opt$input_format, tax_arg)
abund_table  <- ft_obj$abund_table
feature_taxonomy <- ft_obj$feature_taxonomy

# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
message("Loading metadata...")
meta_table <- local({
  sep <- if (grepl("\t", readLines(opt$meta_table, n = 1, warn = FALSE))) "\t" else ","
  read.table(opt$meta_table, header = TRUE, sep = sep, row.names = 1,
             check.names = FALSE, stringsAsFactors = FALSE)
})

# Validate key columns exist before use
.check_col <- function(col, df, arg) {
  if (col != "" && !col %in% colnames(df))
    stop("Column '", col, "' specified by ", arg,
         " not found in metadata. Available columns: ",
         paste(colnames(df), collapse = ", "))
}
.check_col(opt$exclude_column, meta_table, "--exclude_column")
.check_col(opt$type,           meta_table, "--type")
.check_col(opt$type2,          meta_table, "--type2")
if (opt$group != "") {
  for (pc in trimws(strsplit(opt$group, ",")[[1]]))
    .check_col(pc, meta_table, "--group")
}

# ---------------------------------------------------------------------------
# Depth filter
# ---------------------------------------------------------------------------
abund_table <- abund_table[rowSums(abund_table) >= opt$min_library_size, , drop = FALSE]
abund_table <- abund_table[, colSums(abund_table) > 1, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# Align to metadata
abund_table  <- abund_table[rownames(abund_table) %in% rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
meta_table   <- meta_table[rownames(abund_table), , drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# Hypothesis space (exclude, groups, type)
# ---------------------------------------------------------------------------
if (opt$exclude_column != "" && opt$exclude_values != "") {
  exc_vals   <- trimws(strsplit(opt$exclude_values, ",")[[1]])
  keep_rows  <- !meta_table[[opt$exclude_column]] %in% exc_vals
  meta_table <- meta_table[keep_rows, , drop = FALSE]
}

if (opt$group != "") {
  cols <- trimws(strsplit(opt$group, ",")[[1]])
  meta_table$Groups <- if (length(cols) == 1) {
    as.factor(as.character(meta_table[[cols]]))
  } else {
    as.factor(do.call(paste, c(meta_table[, cols, drop = FALSE], sep = " ")))
  }
} else {
  meta_table$Groups <- factor("All")
  message("No --group specified — all samples assigned to group 'All'.")
}

if (opt$type != "") {
  meta_table$Type <- as.factor(as.character(meta_table[[opt$type]]))
} else {
  meta_table$Type <- NULL
}

if (opt$type2 != "") {
  meta_table$Type2 <- as.factor(as.character(meta_table[[opt$type2]]))
} else {
  meta_table$Type2 <- NULL
}

abund_table  <- abund_table[rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# Minimum sample check
# ---------------------------------------------------------------------------
if (nrow(abund_table) < 3)
  stop("At least 3 samples are required after filtering; only ",
       nrow(abund_table), " remain.")

# ---------------------------------------------------------------------------
# Tree
# ---------------------------------------------------------------------------
message("Loading phylogenetic tree...")
feature_tree <- tryCatch(
  read.tree(opt$tree_file),
  error = function(e) stop("Failed to read tree file: ", conditionMessage(e))
)
feature_tree$tip.label <- gsub("'", "", feature_tree$tip.label)

# ---------------------------------------------------------------------------
# Subset to top_n_features by abundance
# ---------------------------------------------------------------------------
n_features_available <- ncol(abund_table)
top_n <- min(opt$top_n_features, n_features_available)
if (top_n < opt$top_n_features)
  message("top_n_features (", opt$top_n_features, ") exceeds available features (",
          n_features_available, ") after filtering — using all ", n_features_available, " features.")
abund_table <- abund_table[, order(colSums(abund_table), decreasing = TRUE),
                           drop = FALSE][, seq_len(top_n), drop = FALSE]

# ---------------------------------------------------------------------------
# Prune tree to match features
# ---------------------------------------------------------------------------
tips_in_data <- feature_tree$tip.label %in% colnames(abund_table)
if (sum(tips_in_data) == 0)
  stop("No tree tips match feature names in the feature table after filtering. ",
       "Check that tree tip labels and feature IDs use the same naming convention.")
feature_tree <- drop.tip(feature_tree,
                     feature_tree$tip.label[!feature_tree$tip.label %in% colnames(abund_table)])

# Align table columns to pruned tree
common_features <- intersect(colnames(abund_table), feature_tree$tip.label)
if (length(common_features) == 0)
  stop("No features remain after aligning feature table to tree tips.")
abund_table <- abund_table[, feature_tree$tip.label, drop = FALSE]
abund_table <- as.matrix(abund_table)

message("Running NRI/NTI on ", nrow(abund_table), " samples x ",
        ncol(abund_table), " features.")

# ---------------------------------------------------------------------------
# SES-MPD (→ NRI) and SES-MNTD (→ NTI)
# ---------------------------------------------------------------------------
second_label <- if (opt$abundance_weighted) "Weighted" else "Unweighted"
cop_dist     <- cophenetic(feature_tree)

if (opt$threads >= 2L) {
  job_mpd  <- parallel::mcparallel(ses.mpd(
    abund_table, cop_dist,
    null.model         = opt$null_model,
    abundance.weighted = opt$abundance_weighted,
    runs               = opt$runs,
    iterations         = opt$iterations
  ))
  job_mntd <- parallel::mcparallel(ses.mntd(
    abund_table, cop_dist,
    null.model         = opt$null_model,
    abundance.weighted = opt$abundance_weighted,
    runs               = opt$runs,
    iterations         = opt$iterations
  ))
  ses_results         <- parallel::mccollect(list(job_mpd, job_mntd))
  abund_table.sesmpd  <- ses_results[[1]]
  abund_table.sesmntd <- ses_results[[2]]
} else {
  abund_table.sesmpd <- ses.mpd(
    abund_table, cop_dist,
    null.model         = opt$null_model,
    abundance.weighted = opt$abundance_weighted,
    runs               = opt$runs,
    iterations         = opt$iterations
  )
  abund_table.sesmntd <- ses.mntd(
    abund_table, cop_dist,
    null.model         = opt$null_model,
    abundance.weighted = opt$abundance_weighted,
    runs               = opt$runs,
    iterations         = opt$iterations
  )
}

# NRI = –1 * MPD z-score;  NTI = –1 * MNTD z-score
nri_vals <- -1 * abund_table.sesmpd$mpd.obs.z
nti_vals <- -1 * abund_table.sesmntd$mntd.obs.z

# ---------------------------------------------------------------------------
# CSV 1: Wide NRI/NTI per sample (+ Groups, optional Type/Type2)
# ---------------------------------------------------------------------------
meta_cols <- c("Groups")
if (!is.null(meta_table$Type))  meta_cols <- c(meta_cols, "Type")
if (!is.null(meta_table$Type2)) meta_cols <- c(meta_cols, "Type2")

data_wide <- data.frame(
  NRI    = nri_vals,
  NTI    = nti_vals,
  meta_table[rownames(abund_table), meta_cols, drop = FALSE],
  check.names = FALSE
)

fname_base <- paste0(second_label, "_", opt$label, "_", opt$null_model)

write.csv(data_wide,
          file = file.path(opt$output_dir,
                           paste0("Environmental_Filtering_", fname_base, ".csv")))
message("Written: Environmental_Filtering_", fname_base, ".csv")

# ---------------------------------------------------------------------------
# CSV 2: Long format with metadata
# ---------------------------------------------------------------------------
meta_sub        <- meta_table[rownames(abund_table), meta_cols, drop = FALSE]
meta_sub$sample <- rownames(meta_sub)

df_nri <- data.frame(sample  = rownames(abund_table),
                     value   = nri_vals,
                     measure = "NRI",
                     stringsAsFactors = FALSE)
df_nti <- data.frame(sample  = rownames(abund_table),
                     value   = nti_vals,
                     measure = "NTI",
                     stringsAsFactors = FALSE)
df_long_base <- rbind(df_nri, df_nti)
df_long      <- merge(df_long_base, meta_sub, by = "sample", all.x = TRUE)

col_order <- c("sample", "value", "measure", "Groups")
if ("Type"  %in% colnames(df_long)) col_order <- c(col_order, "Type")
if ("Type2" %in% colnames(df_long)) col_order <- c(col_order, "Type2")
df_long <- df_long[, col_order, drop = FALSE]

write.csv(df_long,
          file = file.path(opt$output_dir,
                           paste0("EF_long_", fname_base, ".csv")),
          row.names = FALSE)
message("Written: EF_long_", fname_base, ".csv")

# ---------------------------------------------------------------------------
# CSV 3: Pairwise statistical tests
# ---------------------------------------------------------------------------
df_pw_input <- data.frame(
  value   = c(nri_vals, nti_vals),
  Groups  = rep(meta_table[rownames(abund_table), "Groups"], 2),
  measure = rep(c("NRI", "NTI"), each = nrow(abund_table)),
  stringsAsFactors = FALSE
)
df_pw_input <- df_pw_input[complete.cases(df_pw_input$value), ]

all_groups   <- unique(as.character(df_pw_input$Groups))
pairwise_rows <- list()

if (opt$test_method != "none" && length(all_groups) >= 2) {
  group_pairs <- combn(all_groups, 2)
  for (meas in unique(df_pw_input$measure)) {
    df_m <- df_pw_input[df_pw_input$measure == meas, ]
    for (l in seq_len(ncol(group_pairs))) {
      g1  <- group_pairs[1, l]
      g2  <- group_pairs[2, l]
      sub <- df_m[df_m$Groups %in% c(g1, g2), ]
      groups_in_sub <- unique(as.character(sub$Groups))
      if (nrow(sub) < 2 || length(groups_in_sub) < 2) {
        pval <- NA_real_
      } else if (opt$test_method == "anova") {
        pval <- tryCatch({
          p <- summary(aov(value ~ Groups, data = sub))[[1]][["Pr(>F)"]][1]
          if (length(p) == 1L) p else NA_real_
        }, error = function(e) NA_real_)
      } else {
        # kruskal
        pval <- tryCatch(
          kruskal.test(value ~ as.factor(Groups), data = sub)$p.value,
          error = function(e) NA_real_
        )
      }
      pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
        measure     = meas,
        group1      = g1,
        group2      = g2,
        pvalue      = pval,
        test_method = opt$test_method,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(pairwise_rows) > 0) {
  df_pairwise <- do.call(rbind, pairwise_rows)
  if (opt$p_adjust_method != "none") {
    df_pairwise$padj <- p.adjust(df_pairwise$pvalue, method = opt$p_adjust_method)
  } else {
    df_pairwise$padj <- df_pairwise$pvalue
  }
} else {
  df_pairwise <- data.frame(
    measure = character(), group1 = character(), group2 = character(),
    pvalue  = numeric(),   test_method = character(), padj = numeric()
  )
}

write.csv(df_pairwise,
          file = file.path(opt$output_dir,
                           paste0("EF_pairwise_", fname_base, ".csv")),
          row.names = FALSE)
message("Written: EF_pairwise_", fname_base, ".csv")
message("Environmental filtering analysis complete.")
