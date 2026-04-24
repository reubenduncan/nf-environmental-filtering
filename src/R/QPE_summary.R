# QPE_summary.R
# Post-process QPE outputs: assembly process classification and EMS typing.
# Reads the eight CSVs from QPE.R; writes two summary CSVs.
# No phyloseq, no feature table — reads only QPE CSVs.

suppressPackageStartupMessages({
  library(optparse)
})

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--pairwise_bnti_csv",
              type="character", default=NULL,
              help="Path to PairwisebNTI CSV from QPE.R [required]"),
  make_option("--pairwise_rc_csv",
              type="character", default=NULL,
              help="Path to PairwiseRC CSV from QPE.R [required]"),
  make_option("--rc_csv",
              type="character", default=NULL,
              help="Path to RC (summary) CSV from QPE.R [required]"),
  make_option("--coherence_csv",
              type="character", default=NULL,
              help="Path to Coherence CSV from QPE.R [required]"),
  make_option("--boundary_csv",
              type="character", default=NULL,
              help="Path to Boundary CSV from QPE.R [required]"),
  make_option("--turnover_csv",
              type="character", default=NULL,
              help="Path to Turnover CSV from QPE.R [required]"),
  make_option("--output_dir",
              type="character", default=".",
              help="Directory for output CSVs [default: .]"),
  make_option("--label",
              type="character", default="analysis",
              help="Label appended to output filenames [default: analysis]"),
  make_option("--ordering",
              type="character", default="",
              help="Comma-separated group names for ordering output (leave empty for alphabetical)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
required_opts <- c("pairwise_bnti_csv", "pairwise_rc_csv", "rc_csv",
                   "coherence_csv", "boundary_csv", "turnover_csv")
for (r in required_opts) {
  if (is.null(opt[[r]]) || opt[[r]] == "")
    stop("--", r, " is required.")
  if (!file.exists(opt[[r]]))
    stop("File not found for --", r, ": ", opt[[r]])
}

if (!dir.exists(opt$output_dir))
  dir.create(opt$output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
message("Loading QPE result CSVs...")
PairwisebNTI <- read.csv(opt$pairwise_bnti_csv, header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE)
PairwiseRC   <- read.csv(opt$pairwise_rc_csv,   header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE)
RC           <- read.csv(opt$rc_csv,            header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE)
Coherence    <- read.csv(opt$coherence_csv,     header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE)
BoundaryClump <- read.csv(opt$boundary_csv,     header = TRUE, row.names = 1,
                           stringsAsFactors = FALSE)
Turnover     <- read.csv(opt$turnover_csv,      header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Determine group ordering
# ---------------------------------------------------------------------------
if (opt$ordering != "") {
  ordering <- trimws(strsplit(opt$ordering, ",")[[1]])
} else {
  ordering <- sort(unique(as.character(PairwiseRC$Groups)))
}

# ---------------------------------------------------------------------------
# QPE assembly process classification
# Thresholds: |bNTI| > 2 → selection; RC > 0.95 → dispersal limitation;
#             RC < -0.95 → homogenizing dispersal; otherwise undominated.
# ---------------------------------------------------------------------------
message("Classifying assembly processes...")

# Merge pairwise bNTI and RC by row position (both should be same structure)
# Both come from the same loop in QPE.R so rows correspond
if (nrow(PairwisebNTI) != nrow(PairwiseRC))
  warning("PairwisebNTI and PairwiseRC have different row counts (",
          nrow(PairwisebNTI), " vs ", nrow(PairwiseRC),
          "). Results may be unreliable.")

# Build a unified table with complete cases
QPE_table <- data.frame(
  Var1   = as.character(PairwiseRC$Var1),
  Var2   = as.character(PairwiseRC$Var2),
  bNTI   = as.numeric(PairwisebNTI$value),
  RC     = as.numeric(PairwiseRC$value),
  Groups = as.character(PairwiseRC$Groups),
  stringsAsFactors = FALSE
)
QPE_table <- QPE_table[complete.cases(QPE_table), ]
QPE_table$Groups <- factor(QPE_table$Groups)

QPE_df <- NULL
for (grp in levels(QPE_table$Groups)) {
  tmp   <- QPE_table[QPE_table$Groups == grp, ]
  total <- nrow(tmp)
  if (total == 0) next

  sig_selection <- abs(tmp$bNTI) > 2
  hs_count <- sum(tmp$bNTI[sig_selection] < 0)    # bNTI < -2
  vs_count <- sum(tmp$bNTI[sig_selection] > 0)    # bNTI >  2
  nonsig   <- tmp[!sig_selection, ]
  dl_count <- sum(nonsig$RC >  0.95)              # dispersal limitation
  hd_count <- sum(nonsig$RC < -0.95)              # homogenizing dispersal
  ed_count <- nrow(nonsig) - dl_count - hd_count  # undominated

  grp_df <- data.frame(
    variable   = c("Homogeneous Selection", "Variable Selection",
                   "Dispersal Limitation",  "Homogenizing Dispersal",
                   "Undominated"),
    percentage = c(hs_count, vs_count, dl_count, hd_count, ed_count) / total * 100,
    count      = c(hs_count, vs_count, dl_count, hd_count, ed_count),
    total      = total,
    Groups     = grp,
    stringsAsFactors = FALSE
  )
  QPE_df <- if (is.null(QPE_df)) grp_df else rbind(QPE_df, grp_df)
}

if (!is.null(QPE_df)) {
  QPE_df$Groups <- factor(QPE_df$Groups, levels = ordering)
  QPE_df        <- QPE_df[order(QPE_df$Groups), ]
}

write.csv(QPE_df,
          file = file.path(opt$output_dir,
                           paste0("QPE_assembly_processes_", opt$label, ".csv")),
          row.names = FALSE)
message("Written: QPE_assembly_processes_", opt$label, ".csv")

# ---------------------------------------------------------------------------
# EMS classification
# Logic follows Leibold & Mikkelson (2002) and Presley et al. (2010):
#   Coherence z < -1.96 → Checkerboard
#   Coherence z > +1.96 → coherent; then check Turnover:
#     Turnover z < -1.96 → Nested (clumped/hyperdispersed/random species loss)
#     Turnover z > +1.96 → Structured (Clementsian/Gleasonian/Evenly spaced)
#   Coherence p >= 0.05 → Random
# ---------------------------------------------------------------------------
message("Classifying EMS metacommunity types...")

# Helper: detect column by partial name match
.find_col <- function(df, pattern) {
  cols <- grep(pattern, colnames(df), value = TRUE, ignore.case = TRUE)
  if (length(cols) == 0) return(NULL)
  cols[1]
}

# Identify relevant columns (robust to variable column naming from Metacommunity())
coh_z_col <- .find_col(Coherence, "z")
coh_p_col <- .find_col(Coherence, "^p$|pval|Pr")
tur_z_col <- .find_col(Turnover,  "z")
tur_p_col <- .find_col(Turnover,  "^p$|pval|Pr")
bnd_p_col <- .find_col(BoundaryClump, "^p$|pval|Pr")
bnd_i_col <- .find_col(BoundaryClump, "index|clump")

collated_community_types <- character(nrow(Coherence))

for (idx in seq_len(nrow(Coherence))) {
  grp_nm <- rownames(Coherence)[idx]

  # Safe accessors with fallback NA
  coh_z <- if (!is.null(coh_z_col)) as.numeric(Coherence[idx, coh_z_col]) else NA_real_
  coh_p <- if (!is.null(coh_p_col)) as.numeric(Coherence[idx, coh_p_col]) else NA_real_
  tur_z <- if (!is.null(tur_z_col) && grp_nm %in% rownames(Turnover))
              as.numeric(Turnover[grp_nm, tur_z_col]) else NA_real_
  tur_p <- if (!is.null(tur_p_col) && grp_nm %in% rownames(Turnover))
              as.numeric(Turnover[grp_nm, tur_p_col]) else NA_real_
  bnd_p <- if (!is.null(bnd_p_col) && grp_nm %in% rownames(BoundaryClump))
              as.numeric(BoundaryClump[grp_nm, bnd_p_col]) else NA_real_
  bnd_i <- if (!is.null(bnd_i_col) && grp_nm %in% rownames(BoundaryClump))
              as.numeric(BoundaryClump[grp_nm, bnd_i_col]) else NA_real_

  community_type <- "Random"

  if (!is.na(coh_p) && coh_p < 0.05) {
    if (!is.na(coh_z) && coh_z < -1.96) {
      community_type <- "Checkerboard"

    } else if (!is.na(coh_z) && coh_z > 1.96) {
      if (!is.na(tur_z) && tur_z < -1.96) {
        # Nested — boundary clump determines sub-type
        if (!is.na(bnd_p) && bnd_p < 0.05) {
          community_type <- if (!is.na(bnd_i) && bnd_i < 0)
            "Nested Hyperdispersed species loss"
          else
            "Nested Clumped species loss"
        } else {
          community_type <- "Nested Random species loss"
        }
      } else if (!is.na(tur_z) && tur_z > 1.96) {
        # Structured — boundary clump determines sub-type
        if (!is.na(bnd_p) && bnd_p < 0.05) {
          community_type <- if (!is.na(bnd_i) && bnd_i < 0)
            "Evenly spaced"
          else
            "Clementsian"
        } else {
          community_type <- "Gleasonian"
        }
      }
      # Quasi-structure flag if turnover p > 0.05
      if (!is.na(tur_p) && tur_p > 0.05)
        community_type <- paste("Quasi-structure", community_type)
    }
  }
  collated_community_types[idx] <- community_type
}

# Build EMS output — standardize column names regardless of Metacommunity() version
EMS <- data.frame(
  row.names           = rownames(Coherence),
  Coherence           = Coherence,
  Turnover            = Turnover[match(rownames(Coherence), rownames(Turnover)), , drop = FALSE],
  BoundaryClump       = BoundaryClump[match(rownames(Coherence), rownames(BoundaryClump)), , drop = FALSE],
  Metacommunity_type  = collated_community_types,
  stringsAsFactors    = FALSE
)

# Flatten if merge introduced list-columns
EMS <- as.data.frame(lapply(EMS, function(x) if (is.list(x)) unlist(x) else x),
                     stringsAsFactors = FALSE)

write.csv(EMS,
          file = file.path(opt$output_dir, paste0("EMS_", opt$label, ".csv")))
message("Written: EMS_", opt$label, ".csv")
message("QPE summary complete.")
