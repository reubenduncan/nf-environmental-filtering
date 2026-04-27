# environmental-filtering

A Nextflow pipeline for quantifying the relative contributions of stochastic and deterministic processes to microbiome community assembly, using Environmental Filtering (EF), Null model Stochasticity Testing (NST), and Quantitative Process Estimation (QPE).

## Introduction

This pipeline implements three complementary approaches to community assembly analysis:

- **EF (Environmental Filtering):** Mantel-based tests correlating community dissimilarity with environmental distance matrices, identifying which environmental variables drive community structure.
- **NST (Null model Stochasticity Testing):** Compares observed turnover to randomised null distributions to partition stochastic versus deterministic assembly across group pairs.
- **QPE (Quantitative Process Estimation):** Uses phylogenetic null models (βNTI, RC-Bray) to estimate the contributions of homogeneous selection, variable selection, dispersal limitation, homogenising dispersal, and drift.

An optional `--merge_parquet` flag consolidates all 16 output CSVs into a single Parquet file.

## Quick start

```bash
nextflow run main.nf \
  --feature_table /path/to/feature_table.biom \
  --meta_table    /path/to/meta_table.csv \
  --tree_file     /path/to/tree.nwk \
  --groups_column Treatment \
  --label         my_analysis
```

## Parameters

### Input / Output

| Parameter | Default | Description |
|---|---|---|
| `--feature_table` | *(required)* | Path to feature table (BIOM, TSV, or GTDB format) |
| `--meta_table` | *(required)* | Path to sample metadata CSV (first column = sample IDs) |
| `--taxonomy_table` | `""` | Taxonomy TSV (required for `tsv`/`gtdb` input formats) |
| `--tree_file` | `""` | Newick phylogenetic tree (required for QPE/βNTI) |
| `--input_format` | `biom` | `biom` \| `tsv` \| `gtdb` |
| `--output_dir` | `results/` | Directory for output files |

### Filtering

| Parameter | Default | Description |
|---|---|---|
| `--min_library_size` | `5000` | Minimum per-sample read depth; samples below this are dropped |
| `--exclude_column` | `""` | Metadata column used to identify samples for exclusion |
| `--exclude_values` | `""` | Comma-separated values in `exclude_column` to remove |

### Grouping

| Parameter | Default | Description |
|---|---|---|
| `--groups_column` | `""` | Metadata column for the primary grouping variable |
| `--groups_paste_columns` | `""` | Comma-separated columns pasted together to form groups |
| `--type_column` | `""` | Optional secondary metadata column for subsetting |
| `--type2_column` | `""` | Optional tertiary metadata column for subsetting |
| `--type2_levels` | `""` | Comma-separated levels of `type2_column` to retain |

### Environmental Filtering (EF)

| Parameter | Default | Description |
|---|---|---|
| `--ef_variables` | `""` | Comma-separated metadata columns to use as environmental variables |
| `--ef_taxon_rank` | `Feature` | Taxonomic level for EF analysis |
| `--ef_dist_method` | `bray` | Community dissimilarity metric for EF |

### NST

| Parameter | Default | Description |
|---|---|---|
| `--nst_taxon_rank` | `Feature` | Taxonomic level for NST analysis |
| `--nst_dist_method` | `bray` | Dissimilarity metric for NST |
| `--nst_null_model` | `PF` | Null model type: `PF` \| `RF` \| `PCF` \| `RCF` \| `RNDF` |
| `--nst_randomizations` | `999` | Number of null model randomisations |
| `--nst_run_panova` | `true` | Run pairwise ANOVA on stochasticity ratios |

### QPE

| Parameter | Default | Description |
|---|---|---|
| `--qpe_taxon_rank` | `Feature` | Taxonomic level for QPE analysis |
| `--qpe_beta_reps` | `999` | Number of βNTI null model repetitions |
| `--qpe_ems_sims` | `999` | Number of Metacommunity (EMS) simulations — reduce for speed (e.g. `9`) |
| `--qpe_ordering` | `""` | Comma-separated group order for QPE plots |

### Output options

| Parameter | Default | Description |
|---|---|---|
| `--merge_parquet` | `false` | Merge all output CSVs into a single Parquet file |

## Outputs

All files are written to `--output_dir`.

### EF outputs

| File | Description |
|---|---|
| `Environmental_Filtering_{label}.csv` | Wide-format EF results (Mantel r, p-value per variable × group) |
| `EF_long_{label}.csv` | Long-format EF results |
| `EF_pairwise_{label}.csv` | Pairwise group EF comparisons |

### NST outputs

| File | Description |
|---|---|
| `Stochasticity-Ratios_{label}_VALUES.csv` | Stochasticity ratio per group pair |
| `Stochasticity-Ratios_{label}_PAIRWISE.csv` | Pairwise group comparisons of stochasticity |
| `Stochasticity-Ratios_{label}_PANOVA.csv` | Pairwise ANOVA results (produced when `--nst_run_panova true`) |

### QPE outputs

| File | Description |
|---|---|
| `PairwisebNTI_{label}.csv` | Pairwise βNTI values |
| `PairwiseRC_{label}.csv` | Pairwise RC-Bray values |
| `Coherence_{label}.csv` | Metacommunity coherence test results |
| `Boundary_{label}.csv` | Metacommunity boundary clumping results |
| `Turnover_{label}.csv` | Metacommunity species turnover results |
| `Sitescores_{label}.csv` | Site scores from metacommunity ordination |
| `QPE_assembly_processes_{label}.csv` | Estimated contribution (%) of each assembly process per group |
| `EMS_{label}.csv` | EMS z-statistics per group |
| `bNTI_{label}.csv` | βNTI summary statistics per group |
| `RC_{label}.csv` | RC-Bray summary statistics per group |

### Merged output

| File | Description |
|---|---|
| `environmental_filtering_{label}.parquet` | All 16 CSVs merged with `analysis` and `table` metadata columns (`--merge_parquet` only) |

## Performance note

`--qpe_ems_sims` defaults to `999` but each simulation is O(samples × features), making this the most time-consuming step (~45 min/group at default settings on typical datasets). Use `--qpe_ems_sims 9` for exploratory runs.

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/) (default executor — environment built automatically from `environment.yml`)
- **or** Docker with `-profile docker`
- **or** Singularity with `-profile singularity`
- **or** a local R installation with: `optparse`, `vegan`, `ape`, `picante`, `NST`, `metacom`, `stringr`, `data.table`, `phyloseq`, `arrow`

## Running with a local R installation

Add `-profile` to select your execution environment (conda is used by default if no profile is specified):

```bash
nextflow run main.nf \
  -c nextflow.config \
  --feature_table /path/to/table.biom \
  --meta_table    /path/to/meta.csv \
  --tree_file     /path/to/tree.nwk \
  --groups_column Treatment \
  --label         my_analysis \
  --qpe_ems_sims  9
```

Available profiles: `conda` (default), `docker`, `singularity`.