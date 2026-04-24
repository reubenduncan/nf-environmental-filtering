# merge_parquet.R
# Reads all CSVs staged into the Nextflow work directory, stamps each row with
# 'analysis' and 'table' metadata columns, wide-merges with NA-fill for
# differing schemas, and writes a single Parquet file.

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
})

option_list <- list(
  make_option("--label",
              type = "character", default = "analysis",
              help = "Label used in the output filename [default: analysis]"),
  make_option("--output_dir",
              type = "character", default = ".",
              help = "Directory for the output Parquet file [default: .]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$output_dir))
  dir.create(opt$output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Classify each CSV basename into (analysis, table) metadata
# Patterns ordered most-specific-first to avoid prefix collisions.
# ---------------------------------------------------------------------------
classify_csv <- function(fname) {
  b <- basename(fname)
  if      (grepl("^Environmental_Filtering_",              b)) list(analysis = "EF",  table = "ef_wide")
  else if (grepl("^EF_long_",                              b)) list(analysis = "EF",  table = "ef_long")
  else if (grepl("^EF_pairwise_",                          b)) list(analysis = "EF",  table = "ef_pairwise")
  else if (grepl("^Stochasticity-Ratios_.*_VALUES\\.csv$", b)) list(analysis = "NST", table = "nst_values")
  else if (grepl("^Stochasticity-Ratios_.*_PANOVA\\.csv$", b)) list(analysis = "NST", table = "nst_panova")
  else if (grepl("^Stochasticity-Ratios_.*_PAIRWISE\\.csv$",b))list(analysis = "NST", table = "nst_pairwise")
  else if (grepl("^QPE_assembly_processes_",               b)) list(analysis = "QPE", table = "assembly_processes")
  else if (grepl("^EMS_",                                  b)) list(analysis = "QPE", table = "ems")
  else if (grepl("^PairwisebNTI_",                         b)) list(analysis = "QPE", table = "pairwise_bnti")
  else if (grepl("^PairwiseRC_",                           b)) list(analysis = "QPE", table = "pairwise_rc")
  else if (grepl("^Sitescores_",                           b)) list(analysis = "QPE", table = "sitescores")
  else if (grepl("^Coherence_",                            b)) list(analysis = "QPE", table = "coherence")
  else if (grepl("^Boundary_",                             b)) list(analysis = "QPE", table = "boundary")
  else if (grepl("^Turnover_",                             b)) list(analysis = "QPE", table = "turnover")
  else if (grepl("^RC_",                                   b)) list(analysis = "QPE", table = "rc_summary")
  else if (grepl("^bNTI_",                                 b)) list(analysis = "QPE", table = "bnti_summary")
  else                                                          list(analysis = "unknown", table = "unknown")
}

# ---------------------------------------------------------------------------
# Read all CSVs staged into the working directory
# ---------------------------------------------------------------------------
csv_files <- list.files(".", pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
if (length(csv_files) == 0)
  stop("No CSV files found in working directory.")

message("Found ", length(csv_files), " CSV file(s) to merge.")

frames <- lapply(csv_files, function(f) {
  cls <- classify_csv(f)
  df  <- tryCatch(
    read.csv(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) { warning("Failed to read ", basename(f), ": ", conditionMessage(e)); NULL }
  )
  if (is.null(df)) return(NULL)

  # Strip row-index artifacts produced by write.csv(row.names=TRUE):
  # default: first column has empty name; or re-read produces "X" / "X.1" etc.
  if (ncol(df) > 0 && colnames(df)[1] == "") df <- df[, -1, drop = FALSE]
  idx_cols <- grepl("^X$|^X\\.\\d+$", colnames(df))
  if (any(idx_cols)) df <- df[, !idx_cols, drop = FALSE]

  df$analysis    <- cls$analysis
  df$table       <- cls$table
  df$source_file <- basename(f)

  message("  ", basename(f), " → analysis=", cls$analysis,
          ", table=", cls$table, " (", nrow(df), " rows)")
  df
})

frames <- Filter(Negate(is.null), frames)
if (length(frames) == 0)
  stop("All CSV files failed to read — cannot write Parquet.")

# ---------------------------------------------------------------------------
# Wide merge: collect superset of columns, NA-fill missing, rbind
# ---------------------------------------------------------------------------
all_cols       <- unique(unlist(lapply(frames, colnames)))
frames_aligned <- lapply(frames, function(df) {
  df[setdiff(all_cols, colnames(df))] <- NA
  df[, all_cols, drop = FALSE]
})

merged <- do.call(rbind, frames_aligned)
rownames(merged) <- NULL

# ---------------------------------------------------------------------------
# Write Parquet
# ---------------------------------------------------------------------------
out_path <- file.path(opt$output_dir,
                      paste0("environmental_filtering_", opt$label, ".parquet"))
arrow::write_parquet(merged, out_path)
message("Written: ", out_path,
        " (", nrow(merged), " rows x ", ncol(merged), " columns)")
