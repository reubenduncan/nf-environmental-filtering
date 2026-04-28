# QPE.R
# Quantitative Process Estimation: betaNTI, Raup-Crick (abundance + incidence),
# and Elements of Metacommunity Structure (EMS) per group.
# Outputs eight CSVs consumed by QPE_summary.R.

local({
  conda_prefix <- Sys.getenv("CONDA_PREFIX")
  lib <- if (nchar(conda_prefix) > 0)
    file.path(conda_prefix, "lib", "R", "library")
  else
    .libPaths()[1]
  pkgs <- c("ecodist", "metacom")
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) {
    message("Installing missing R packages into ", lib, ": ",
            paste(miss, collapse = ", "))
    install.packages(miss, lib = lib, repos = "https://cloud.r-project.org",
                     quiet = TRUE)
  }
})

suppressPackageStartupMessages({
  library(optparse)
  library(phyloseq)
  library(ape)
  library(picante)
  library(ecodist)
  library(metacom)
  library(reshape2)
  library(plyr)
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

  # QPE parameters
  make_option("--beta_reps",
              type="integer", default=999L,
              help="Randomizations for betaNTI and Raup-Crick [default: 999]"),
  make_option("--ems_sims",
              type="integer", default=999L,
              help="Null model simulations for metacom EMS [default: 999]")
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
  stop("--tree_file is required for QPE (betaNTI calculation).")
if (!opt$input_format %in% c("biom", "tsv", "gtdb"))
  stop("--input_format must be one of: biom, tsv, gtdb.")

if (!dir.exists(opt$output_dir))
  dir.create(opt$output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Helper: Summarize a distance matrix column (mean, sd, N, se, ci) by group
# ---------------------------------------------------------------------------
summarySE <- function(data, measurevar, groupvars, na.rm = TRUE, conf.interval = 0.95) {
  datac <- ddply(data, groupvars, .drop = TRUE,
                 .fun = function(xx, col) {
                   vals <- xx[[col]]
                   if (na.rm) vals <- vals[!is.na(vals)]
                   n    <- length(vals)
                   mn   <- if (n > 0) mean(vals) else NA_real_
                   s    <- if (n > 1) sd(vals)   else NA_real_
                   c(N = n, mean = mn, sd = s)
                 },
                 measurevar)
  colnames(datac)[colnames(datac) == "mean"] <- measurevar
  datac$se <- datac$sd / sqrt(datac$N)
  datac$ci <- datac$se * qt(conf.interval / 2 + 0.5, datac$N - 1)
  datac
}

# ---------------------------------------------------------------------------
# Raup-Crick (incidence-based)
# ---------------------------------------------------------------------------
raup_crick <- function(spXsite,
                       plot_names_in_col1   = FALSE,
                       classic_metric       = FALSE,
                       split_ties           = TRUE,
                       reps                 = 999,
                       set_all_species_equal = FALSE,
                       as.distance.matrix   = TRUE,
                       report_similarity    = FALSE) {

  if (plot_names_in_col1) {
    row.names(spXsite) <- spXsite[, 1]
    spXsite <- spXsite[, -1]
  }

  n_sites <- nrow(spXsite)
  gamma   <- ncol(spXsite)

  # Convert to incidence
  spXsite_inc <- ceiling(spXsite / max(spXsite))

  occur <- apply(spXsite_inc, MARGIN = 2, FUN = sum)
  if (set_all_species_equal) occur <- rep(1, gamma)

  alpha_levels <- sort(unique(apply(spXsite_inc, MARGIN = 1, FUN = sum)))
  alpha_table  <- data.frame(smaller_alpha = NA_real_, bigger_alpha = NA_real_)
  null_array   <- list()
  col_count    <- 1

  for (a1 in seq_along(alpha_levels)) {
    for (a2 in a1:length(alpha_levels)) {
      null_shared_spp <- numeric(reps)
      for (i in seq_len(reps)) {
        com1 <- rep(0L, gamma)
        com2 <- rep(0L, gamma)
        com1[sample.int(gamma, alpha_levels[a1], replace = FALSE, prob = occur)] <- 1L
        com2[sample.int(gamma, alpha_levels[a2], replace = FALSE, prob = occur)] <- 1L
        null_shared_spp[i] <- sum((com1 + com2) > 1)
      }
      null_array[[col_count]] <- null_shared_spp
      alpha_table[col_count, "smaller_alpha"] <- alpha_levels[a1]
      alpha_table[col_count, "bigger_alpha"]  <- alpha_levels[a2]
      col_count <- col_count + 1
    }
  }
  alpha_table$matching <- paste(alpha_table[, 1], alpha_table[, 2], sep = "_")

  results <- matrix(NA_real_, nrow = n_sites, ncol = n_sites,
                    dimnames = list(rownames(spXsite), rownames(spXsite)))

  for (i in seq_len(n_sites)) {
    for (j in seq_len(n_sites)) {
      if (i == j) next
      n_shared_obs <- sum((spXsite_inc[i, ] + spXsite_inc[j, ]) > 1)
      obs_a_pair   <- sort(c(sum(spXsite_inc[i, ]), sum(spXsite_inc[j, ])))
      null_index   <- which(alpha_table$matching ==
                              paste(obs_a_pair[1], obs_a_pair[2], sep = "_"))
      if (length(null_index) == 0) { results[i, j] <- NA_real_; next }
      nv          <- null_array[[null_index]]
      num_exact   <- sum(nv == n_shared_obs)
      num_greater <- sum(nv >  n_shared_obs)
      rc <- if (split_ties) (num_greater + num_exact / 2) / reps else num_greater / reps
      if (!classic_metric)  rc <- (rc - 0.5) * 2
      if (report_similarity && !classic_metric) rc <- rc * -1
      if (report_similarity &&  classic_metric) rc <- 1 - rc
      results[i, j] <- round(rc, 2)
    }
  }
  if (as.distance.matrix) results <- as.dist(results)
  results
}

# ---------------------------------------------------------------------------
# Raup-Crick abundance-weighted (Bray-Curtis null)
# ---------------------------------------------------------------------------
raup_crick_abundance <- function(spXsite,
                                  plot_names_in_col1    = FALSE,
                                  classic_metric        = FALSE,
                                  split_ties            = TRUE,
                                  reps                  = 9999,
                                  set_all_species_equal = FALSE,
                                  as.distance.matrix    = TRUE,
                                  report_similarity     = FALSE) {

  if (plot_names_in_col1) {
    row.names(spXsite) <- spXsite[, 1]
    spXsite <- spXsite[, -1]
  }

  n_sites     <- nrow(spXsite)
  gamma       <- ncol(spXsite)
  spXsite_inc <- ceiling(spXsite / max(spXsite))
  occur       <- apply(spXsite_inc, MARGIN = 2, FUN = sum)
  abundance   <- apply(spXsite,     MARGIN = 2, FUN = sum)

  results <- matrix(NA_real_, nrow = n_sites, ncol = n_sites,
                    dimnames = list(rownames(spXsite), rownames(spXsite)))

  for (null.one in seq_len(n_sites - 1)) {
    for (null.two in (null.one + 1):n_sites) {
      null_bray_curtis <- numeric(reps)

      alpha_one <- sum(spXsite_inc[null.one, ])
      alpha_two <- sum(spXsite_inc[null.two, ])

      for (i in seq_len(reps)) {
        com1 <- rep(0L, gamma)
        com2 <- rep(0L, gamma)

        # Draw incidence from regional pool weighted by occurrence
        com1[sample.int(gamma, alpha_one, replace = FALSE, prob = occur)] <- 1L
        # Distribute remaining reads proportionally to abundance among selected spp
        extra1 <- sum(spXsite[null.one, ]) - alpha_one
        if (extra1 > 0 && any(com1 > 0)) {
          samp1   <- sample(which(com1 > 0), extra1, replace = TRUE,
                            prob = abundance[which(com1 > 0)])
          counts1 <- tabulate(samp1, nbins = gamma)
          com1    <- com1 + counts1
        }

        com2[sample.int(gamma, alpha_two, replace = FALSE, prob = occur)] <- 1L
        extra2 <- sum(spXsite[null.two, ]) - alpha_two
        if (extra2 > 0 && any(com2 > 0)) {
          samp2   <- sample(which(com2 > 0), extra2, replace = TRUE,
                            prob = abundance[which(com2 > 0)])
          counts2 <- tabulate(samp2, nbins = gamma)
          com2    <- com2 + counts2
        }

        null.spXsite       <- rbind(com1, com2)
        null_bray_curtis[i] <- ecodist::distance(null.spXsite, method = "bray-curtis")
      }

      obs.bray  <- ecodist::distance(spXsite[c(null.one, null.two), ], method = "bray-curtis")
      num_exact <- sum(null_bray_curtis == obs.bray)
      num_less  <- sum(null_bray_curtis <  obs.bray)

      rc <- if (split_ties) (num_less + num_exact / 2) / reps else num_less / reps
      if (!classic_metric) rc <- (rc - 0.5) * 2
      results[null.two, null.one] <- round(rc, 2)

      message("RC_abundance: site pair [", null.one, ",", null.two, "] — ", date())
    }
  }
  if (as.distance.matrix) results <- as.dist(results)
  results
}

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

.check_col <- function(col, df, arg) {
  if (col != "" && !col %in% colnames(df))
    stop("Column '", col, "' specified by ", arg,
         " not found in metadata. Available columns: ",
         paste(colnames(df), collapse = ", "))
}
.check_col(opt$exclude_column,  meta_table, "--exclude_column")
if (opt$group != "") {
  for (pc in trimws(strsplit(opt$group, ",")[[1]]))
    .check_col(pc, meta_table, "--group")
}

# ---------------------------------------------------------------------------
# Depth filter + alignment
# ---------------------------------------------------------------------------
abund_table <- abund_table[rowSums(abund_table) >= opt$min_library_size, , drop = FALSE]
abund_table <- abund_table[, colSums(abund_table) > 1, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

abund_table  <- abund_table[rownames(abund_table) %in% rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
meta_table   <- meta_table[rownames(abund_table), , drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# Hypothesis space
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

abund_table  <- abund_table[rownames(meta_table), , drop = FALSE]
abund_table  <- abund_table[, colSums(abund_table) > 0, drop = FALSE]
feature_taxonomy <- feature_taxonomy[colnames(abund_table), , drop = FALSE]

# ---------------------------------------------------------------------------
# Phylogenetic tree
# ---------------------------------------------------------------------------
message("Loading phylogenetic tree...")
feature_tree <- tryCatch(
  read.tree(opt$tree_file),
  error = function(e) stop("Failed to read tree file: ", conditionMessage(e))
)
feature_tree$tip.label <- gsub("'", "", feature_tree$tip.label)

tips_in_data <- feature_tree$tip.label %in% colnames(abund_table)
if (sum(tips_in_data) == 0)
  stop("No tree tips match feature names in the feature table after filtering.")

# ---------------------------------------------------------------------------
# Build phyloseq and rarefy
# ---------------------------------------------------------------------------
OTU_ps  <- otu_table(as.matrix(abund_table),  taxa_are_rows = FALSE)
TAX_ps  <- tax_table(as.matrix(feature_taxonomy))
SAM_ps  <- sample_data(meta_table)
physeq  <- merge_phyloseq(phyloseq(OTU_ps, TAX_ps), SAM_ps, feature_tree)
physeq  <- prune_taxa(taxa_sums(physeq) > 10, physeq)

min_depth <- min(sample_sums(physeq))
message("Rarefying to depth: ", min_depth)
set.seed(42)
physeq_rel <- rarefy_even_depth(physeq, sample.size = min_depth, verbose = FALSE)

abund_table_ems <- as.matrix(otu_table(physeq_rel))    # samples x features
meta_table_ems  <- as.data.frame(sample_data(physeq_rel))
beta.reps       <- opt$beta_reps

# Refresh tree after rarefaction pruning
feature_tree <- phy_tree(physeq_rel)

# ---------------------------------------------------------------------------
# Collation containers
# ---------------------------------------------------------------------------
collated_coherence              <- NULL
collated_turnover               <- NULL
collated_boundary               <- NULL
collated_sitescores             <- NULL
collated_pairwise_RC_abundance  <- NULL
collated_pairwise_RC_incidence  <- NULL
collated_pairwise_bNTI         <- NULL

grp_levels <- levels(meta_table_ems$Groups)
message("Running QPE analysis over ", length(grp_levels), " group(s): ",
        paste(grp_levels, collapse = ", "))

# ---------------------------------------------------------------------------
# Per-group loop
# ---------------------------------------------------------------------------
for (i in seq_along(grp_levels)) {
  grp_name <- grp_levels[i]
  message("\n=== Group: ", grp_name, " ===")

  grp_mask              <- meta_table_ems$Groups == grp_name
  abund_table_ems_group <- abund_table_ems[grp_mask, , drop = FALSE]
  abund_table_ems_group <- abund_table_ems_group[,
                             colSums(abund_table_ems_group) > 0, drop = FALSE]
  n_samp <- nrow(abund_table_ems_group)
  message("  Samples in group: ", n_samp)

  # ---- EMS (Metacommunity) ------------------------------------------------
  coherence  <- NULL
  turnover   <- NULL
  boundary   <- NULL
  sitescores <- NULL

  tryCatch({
    met_ems <- Metacommunity(abund_table_ems_group,
                             scores      = 1,
                             method      = "r1",
                             sims        = opt$ems_sims,
                             order       = TRUE,
                             binary      = FALSE,
                             verbose     = TRUE,
                             allowEmpty  = TRUE)
    om_ems  <- OrderMatrix(abund_table_ems_group, outputScores = TRUE, binary = FALSE)

    coh_rows <- met_ems$Coherence[-nrow(met_ems$Coherence), , drop = FALSE]
    coherence  <- as.data.frame(t(setNames(coh_rows[, 2], coh_rows[, 1])),
                                stringsAsFactors = FALSE)
    rownames(coherence) <- grp_name

    tur_rows   <- met_ems$Turnover[-nrow(met_ems$Turnover), , drop = FALSE]
    turnover   <- as.data.frame(t(setNames(tur_rows[, 2], tur_rows[, 1])),
                                stringsAsFactors = FALSE)
    rownames(turnover) <- grp_name

    bnd_rows   <- met_ems$Boundary
    boundary   <- as.data.frame(t(setNames(bnd_rows[, 2], bnd_rows[, 1])),
                                stringsAsFactors = FALSE)
    rownames(boundary) <- grp_name

    sitescores <- data.frame(sitescores = om_ems$sitescores, Groups = grp_name,
                             stringsAsFactors = FALSE)
  }, error = function(e) {
    warning("Metacommunity() failed for group '", grp_name, "': ", conditionMessage(e))
  })

  # ---- Raup-Crick (abundance) ---------------------------------------------
  pairwise_RC_abundance <- NULL
  tryCatch({
    rc_abund <- raup_crick_abundance(abund_table_ems_group,
                                     set_all_species_equal = FALSE,
                                     plot_names_in_col1    = FALSE,
                                     reps                  = beta.reps)
    pairwise_RC_abundance <- reshape2::melt(as.matrix(rc_abund))
    pairwise_RC_abundance$Groups <- grp_name
  }, error = function(e) {
    warning("raup_crick_abundance() failed for group '", grp_name, "': ", conditionMessage(e))
  })

  # ---- Raup-Crick (incidence) ---------------------------------------------
  pairwise_RC_incidence <- NULL
  tryCatch({
    rc_incid <- raup_crick(abund_table_ems_group,
                            plot_names_in_col1    = FALSE,
                            reps                  = beta.reps,
                            as.distance.matrix    = TRUE,
                            set_all_species_equal = FALSE)
    pairwise_RC_incidence <- reshape2::melt(as.matrix(rc_incid))
    pairwise_RC_incidence$Groups <- grp_name
  }, error = function(e) {
    warning("raup_crick() failed for group '", grp_name, "': ", conditionMessage(e))
  })

  # ---- betaNTI ------------------------------------------------------------
  pairwise_bNTI <- NULL
  if (n_samp < 3) {
    warning("Group '", grp_name, "' has only ", n_samp, " sample(s); ",
            "skipping betaNTI (requires >= 3 samples for meaningful pairwise comparisons).")
  } else {
    tryCatch({
      m_phylo    <- match.phylo.data(feature_tree, t(abund_table_ems_group))
      tree_grp   <- m_phylo$phy
      at_grp     <- t(abund_table_ems_group)   # features x samples

      cop_grp    <- cophenetic(tree_grp)
      beta_obs   <- as.matrix(comdistnt(t(at_grp), cop_grp, abundance.weighted = TRUE))

      rand_bMNTD <- array(-999, dim = c(ncol(at_grp), ncol(at_grp), beta.reps))
      for (rep in seq_len(beta.reps)) {
        rand_bMNTD[, , rep] <- as.matrix(
          comdistnt(t(at_grp),
                    taxaShuffle(cop_grp),
                    abundance.weighted    = TRUE,
                    exclude.conspecifics = FALSE))
        message("  betaNTI rep ", rep, "/", beta.reps, " — ", date())
      }

      weighted_bNTI <- matrix(NA_real_, nrow = ncol(at_grp), ncol = ncol(at_grp),
                              dimnames = list(colnames(at_grp), colnames(at_grp)))
      for (cols in seq_len(ncol(at_grp) - 1)) {
        for (rows in (cols + 1):ncol(at_grp)) {
          rv <- rand_bMNTD[rows, cols, ]
          rv <- rv[rv != -999]
          if (length(rv) < 2) next
          weighted_bNTI[rows, cols] <- (beta_obs[rows, cols] - mean(rv)) / sd(rv)
        }
      }
      pairwise_bNTI <- reshape2::melt(as.matrix(weighted_bNTI))
      pairwise_bNTI$Groups <- grp_name
    }, error = function(e) {
      warning("betaNTI failed for group '", grp_name, "': ", conditionMessage(e))
    })
  }

  # ---- Collate ------------------------------------------------------------
  .append <- function(acc, x) if (is.null(acc)) x else rbind(acc, x)
  if (!is.null(coherence))             collated_coherence             <- .append(collated_coherence,             coherence)
  if (!is.null(turnover))              collated_turnover              <- .append(collated_turnover,              turnover)
  if (!is.null(boundary))             collated_boundary              <- .append(collated_boundary,              boundary)
  if (!is.null(sitescores))            collated_sitescores            <- .append(collated_sitescores,            sitescores)
  if (!is.null(pairwise_RC_abundance)) collated_pairwise_RC_abundance <- .append(collated_pairwise_RC_abundance, pairwise_RC_abundance)
  if (!is.null(pairwise_RC_incidence)) collated_pairwise_RC_incidence <- .append(collated_pairwise_RC_incidence, pairwise_RC_incidence)
  if (!is.null(pairwise_bNTI))         collated_pairwise_bNTI         <- .append(collated_pairwise_bNTI,         pairwise_bNTI)
}

# ---------------------------------------------------------------------------
# Summarize RC and bNTI across groups
# ---------------------------------------------------------------------------
collated_RC <- if (!is.null(collated_pairwise_RC_incidence)) summarySE(collated_pairwise_RC_incidence, measurevar = "value", groupvars = "Groups") else NULL

collated_bNTI <- if (!is.null(collated_pairwise_bNTI)) summarySE(collated_pairwise_bNTI, measurevar = "value", groupvars = "Groups") else NULL

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
lbl <- opt$label

.safe_write <- function(obj, path) {
  if (!is.null(obj)) {
    write.csv(obj, file = path)
    message("Written: ", basename(path))
  } else {
    message("Skipped (no data): ", basename(path))
  }
}

.safe_write(collated_coherence,
            file.path(opt$output_dir, paste0("Coherence_",    lbl, ".csv")))
.safe_write(collated_boundary,
            file.path(opt$output_dir, paste0("Boundary_",     lbl, ".csv")))
.safe_write(collated_turnover,
            file.path(opt$output_dir, paste0("Turnover_",     lbl, ".csv")))
.safe_write(collated_sitescores,
            file.path(opt$output_dir, paste0("Sitescores_",   lbl, ".csv")))
.safe_write(collated_pairwise_RC_abundance,
            file.path(opt$output_dir, paste0("PairwiseRC_",   lbl, ".csv")))
.safe_write(collated_pairwise_bNTI,
            file.path(opt$output_dir, paste0("PairwisebNTI_", lbl, ".csv")))

if (!is.null(collated_RC)) {
  rownames(collated_RC) <- collated_RC[, 1]
  collated_RC           <- collated_RC[, -1, drop = FALSE]
  .safe_write(collated_RC,
              file.path(opt$output_dir, paste0("RC_",   lbl, ".csv")))
}
if (!is.null(collated_bNTI)) {
  rownames(collated_bNTI) <- collated_bNTI[, 1]
  collated_bNTI           <- collated_bNTI[, -1, drop = FALSE]
  .safe_write(collated_bNTI,
              file.path(opt$output_dir, paste0("bNTI_", lbl, ".csv")))
}

message("QPE analysis complete.")
