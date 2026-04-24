# NST.R
# Null Stochasticity Test (NST/ST/MST) for microbiome assembly processes.
# Outputs three CSVs: values per group, PANOVA results, pairwise values.

suppressPackageStartupMessages({
  library(optparse)
  library(NST)
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
  make_option("--groups_column",
              type="character", default="",
              help="Metadata column for the Groups variable"),
  make_option("--groups_paste_columns",
              type="character", default="",
              help="Comma-separated metadata columns pasted together to form Groups"),

  # NST parameters
  make_option("--nst_randomizations",
              type="integer", default=1000L,
              help="Number of null model randomizations [default: 1000]"),
  make_option("--nst_distance",
              type="character", default="cao",
              help="Distance measure: bray | jaccard | cao | chao [default: cao]"),
  make_option("--nst_null_model",
              type="character", default="PF",
              help="Null model: EE|EP|EF|PE|PP|PF|FE|FP|FF [default: PF]"),
  make_option("--nst_abundance_weighted",
              action="store_true", default=FALSE,
              help="Use abundance-weighted (Ruzicka) instead of incidence-based (Jaccard) [default: FALSE]"),
  make_option("--nst_ses",
              action="store_true", default=FALSE,
              help="Calculate standardized effect size [default: TRUE via nextflow.config]"),
  make_option("--nst_rc",
              action="store_true", default=FALSE,
              help="Calculate modified Raup-Crick metric [default: FALSE]")
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
if (!opt$input_format %in% c("biom", "tsv", "gtdb"))
  stop("--input_format must be one of: biom, tsv, gtdb.")
if (!opt$nst_distance %in% c("bray", "jaccard", "cao", "chao"))
  stop("--nst_distance must be one of: bray, jaccard, cao, chao.")
valid_null_models <- c("EE","EP","EF","PE","PP","PF","FE","FP","FF")
if (!opt$nst_null_model %in% valid_null_models)
  stop("--nst_null_model must be one of: ", paste(valid_null_models, collapse=", "))

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
OTU_taxonomy <- ft_obj$OTU_taxonomy

# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
message("Loading metadata...")
meta_table <- read.csv(opt$meta_table, header = TRUE, row.names = 1,
                       stringsAsFactors = FALSE)

.check_col <- function(col, df, arg) {
  if (col != "" && !col %in% colnames(df))
    stop("Column '", col, "' specified by ", arg,
         " not found in metadata. Available columns: ",
         paste(colnames(df), collapse = ", "))
}
.check_col(opt$exclude_column,  meta_table, "--exclude_column")
.check_col(opt$groups_column,   meta_table, "--groups_column")
if (opt$groups_paste_columns != "") {
  for (pc in trimws(strsplit(opt$groups_paste_columns, ",")[[1]]))
    .check_col(pc, meta_table, "--groups_paste_columns")
}

# ---------------------------------------------------------------------------
# Depth filter + alignment
# ---------------------------------------------------------------------------
abund_table <- abund_table[rowSums(abund_table) >= opt$min_library_size, , drop = FALSE]
abund_table <- abund_table[, colSums(abund_table) > 1, drop = FALSE]
OTU_taxonomy <- OTU_taxonomy[colnames(abund_table), , drop = FALSE]

abund_table  <- abund_table[rownames(abund_table) %in% rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
meta_table   <- meta_table[rownames(abund_table), , drop = FALSE]
OTU_taxonomy <- OTU_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# Hypothesis space
# ---------------------------------------------------------------------------
if (opt$exclude_column != "" && opt$exclude_values != "") {
  exc_vals   <- trimws(strsplit(opt$exclude_values, ",")[[1]])
  keep_rows  <- !meta_table[[opt$exclude_column]] %in% exc_vals
  meta_table <- meta_table[keep_rows, , drop = FALSE]
}

if (opt$groups_paste_columns != "") {
  paste_cols        <- trimws(strsplit(opt$groups_paste_columns, ",")[[1]])
  meta_table$Groups <- as.factor(do.call(paste, c(meta_table[, paste_cols, drop = FALSE], sep = " ")))
} else if (opt$groups_column != "") {
  meta_table$Groups <- as.factor(as.character(meta_table[[opt$groups_column]]))
} else {
  meta_table$Groups <- factor("All")
}

abund_table  <- abund_table[rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
OTU_taxonomy <- OTU_taxonomy[colnames(abund_table), , drop = FALSE]
abund_mat    <- as(abund_table, "matrix")

# ---------------------------------------------------------------------------
# Validate group structure: need >= 2 groups each with >= 2 samples
# ---------------------------------------------------------------------------
grp_counts <- table(meta_table$Groups)
valid_grps <- names(grp_counts)[grp_counts >= 2]
if (length(valid_grps) < 2)
  stop("NST requires at least 2 groups with >= 2 samples each. ",
       "Found ", length(valid_grps), " qualifying group(s) after filtering. ",
       "Check --groups_column, --min_library_size, and --exclude_values.")

if (length(valid_grps) < nlevels(meta_table$Groups)) {
  dropped <- setdiff(levels(meta_table$Groups), valid_grps)
  warning("Groups with < 2 samples dropped from NST analysis: ",
          paste(dropped, collapse = ", "))
  keep_samp     <- meta_table$Groups %in% valid_grps
  meta_table    <- meta_table[keep_samp, , drop = FALSE]
  abund_mat     <- abund_mat[keep_samp, , drop = FALSE]
  meta_table$Groups <- droplevels(meta_table$Groups)
}

message("Running tNST on ", nrow(abund_mat), " samples across ",
        nlevels(meta_table$Groups), " groups.")

# ---------------------------------------------------------------------------
# NST analysis
# ---------------------------------------------------------------------------
tnst <- tryCatch(
  tNST(comm        = abund_mat,
       group       = meta_table[, "Groups", drop = FALSE],
       rand        = opt$nst_randomizations,
       dist.method = opt$nst_distance,
       null.model  = opt$nst_null_model,
       output.rand = TRUE,
       nworker     = 1,
       SES         = opt$nst_ses,
       RC          = opt$nst_rc),
  error = function(e) {
    stop("tNST() failed: ", conditionMessage(e))
  }
)

# ---------------------------------------------------------------------------
# Build per-group long data frame (NST / ST / MST)
# ---------------------------------------------------------------------------
df_values <- NULL
for (metric_prefix in c("NST", "ST", "MST")) {
  col_match <- names(tnst$index.grp)[grepl(paste0("^", metric_prefix), names(tnst$index.grp))]
  if (length(col_match) == 0) next
  tmp <- tnst$index.grp[, c("group", col_match[1]), drop = FALSE]
  colnames(tmp) <- c("Groups", "value")
  tmp$measure   <- metric_prefix
  df_values <- if (is.null(df_values)) tmp else rbind(df_values, tmp)
}

fname_base <- paste(opt$nst_null_model, opt$nst_distance,
                    opt$nst_randomizations, opt$nst_ses, opt$nst_rc,
                    opt$nst_abundance_weighted, opt$label, sep = "_")

# ---------------------------------------------------------------------------
# CSV 1: NST/ST/MST values per group (long format)
# ---------------------------------------------------------------------------
write.csv(df_values,
          file = file.path(opt$output_dir,
                           paste0("Stochasticity-Ratios_", fname_base, "_VALUES.csv")),
          row.names = FALSE)
message("Written: Stochasticity-Ratios_", fname_base, "_VALUES.csv")

# ---------------------------------------------------------------------------
# CSV 2: PANOVA results
# ---------------------------------------------------------------------------
nst_pova <- tryCatch(
  nst.panova(nst.result = tnst, rand = opt$nst_randomizations),
  error = function(e) {
    warning("nst.panova() failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(nst_pova)) {
  # Replace integer group codes with group names
  grp_levels <- levels(meta_table$Groups)
  for (col_nm in c("group1", "group2")) {
    if (col_nm %in% names(nst_pova)) {
      for (i in seq_along(grp_levels))
        nst_pova[[col_nm]] <- gsub(paste0("^", i, "$"), grp_levels[i], nst_pova[[col_nm]])
    }
  }
  write.csv(nst_pova,
            file = file.path(opt$output_dir,
                             paste0("Stochasticity-Ratios_", fname_base, "_PANOVA.csv")))
  message("Written: Stochasticity-Ratios_", fname_base, "_PANOVA.csv")
} else {
  message("PANOVA not written due to error in nst.panova().")
}

# ---------------------------------------------------------------------------
# CSV 3: Pairwise values
# ---------------------------------------------------------------------------
write.csv(tnst$index.pair.grp,
          file = file.path(opt$output_dir,
                           paste0("Stochasticity-Ratios_", fname_base, "_PAIRWISE.csv")))
message("Written: Stochasticity-Ratios_", fname_base, "_PAIRWISE.csv")
message("NST analysis complete.")
