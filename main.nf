#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Environmental Filtering Pipeline
// Runs NRI/NTI, NST stochasticity ratios, and QPE/metacommunity analyses
// in parallel, then post-processes QPE outputs.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Process: ENV_FILTERING — NRI/NTI via picante
// ---------------------------------------------------------------------------
process ENV_FILTERING {
    tag "${params.label}"

    publishDir "${params.output_dir}/EF", mode: 'copy'

    input:
    path feature_table
    path meta_table
    path tree_file

    output:
    path "Environmental_Filtering_*.csv", emit: ef_wide
    path "EF_long_*.csv",                emit: ef_long
    path "EF_pairwise_*.csv",            emit: ef_pairwise

    script:
    def tax_arg  = params.taxonomy_table       ? "--taxonomy_table '${params.taxonomy_table}'" : ""
    def excl_col = params.exclude_column       ? "--exclude_column '${params.exclude_column}'" : ""
    def excl_val = params.exclude_values       ? "--exclude_values '${params.exclude_values}'" : ""
    def aw_arg   = params.ef_abundance_weighted ? "--abundance_weighted"                        : ""

    """
    Rscript ${projectDir}/src/R/environmental_filtering.R \\
        --feature_table     '${feature_table}'             \\
        --meta_table        '${meta_table}'                \\
        --tree_file         '${tree_file}'                 \\
        --input_format      '${params.input_format}'       \\
        --output_dir        '.'                            \\
        --label             '${params.label}'              \\
        --scripts_dir       '${projectDir}'                \\
        --min_library_size  ${params.min_library_size}     \\
        --runs              ${params.ef_runs}              \\
        --iterations        ${params.ef_iterations}        \\
        --top_n_features    ${params.ef_top_n_features}    \\
        --null_model        '${params.ef_null_model}'      \\
        --test_method       '${params.ef_test_method}'     \\
        --p_adjust_method   '${params.ef_p_adjust_method}' \\
        ${aw_arg}                                          \\
        ${tax_arg}                                         \\
        ${excl_col}                                        \\
        ${excl_val}                                        \\
        ${params.group ? "--group '${params.group}'" : ""} \\
        ${params.type  ? "--type  '${params.type}'"  : ""} \\
        ${params.type2 ? "--type2 '${params.type2}'" : ""}
    """
}

// ---------------------------------------------------------------------------
// Process: NST — Null stochasticity test
// ---------------------------------------------------------------------------
process NST {
    tag "${params.label}"

    publishDir "${params.output_dir}/NST", mode: 'copy'

    input:
    path feature_table
    path meta_table

    output:
    path "Stochasticity-Ratios_*_VALUES.csv",   emit: nst_values
    path "Stochasticity-Ratios_*_PANOVA.csv",   emit: nst_panova,   optional: true
    path "Stochasticity-Ratios_*_PAIRWISE.csv", emit: nst_pairwise

    script:
    def tax_arg  = params.taxonomy_table        ? "--taxonomy_table '${params.taxonomy_table}'" : ""
    def excl_col = params.exclude_column        ? "--exclude_column '${params.exclude_column}'" : ""
    def excl_val = params.exclude_values        ? "--exclude_values '${params.exclude_values}'" : ""
    def aw_arg   = params.nst_abundance_weighted ? "--nst_abundance_weighted"                   : ""
    def ses_arg  = params.nst_ses               ? "--nst_ses"                                   : ""
    def rc_arg   = params.nst_rc                ? "--nst_rc"                                    : ""

    """
    Rscript ${projectDir}/src/R/NST.R \\
        --feature_table        '${feature_table}'            \\
        --meta_table           '${meta_table}'               \\
        --input_format         '${params.input_format}'      \\
        --output_dir           '.'                           \\
        --label                '${params.label}'             \\
        --scripts_dir          '${projectDir}'               \\
        --min_library_size     ${params.min_library_size}    \\
        --nst_randomizations   ${params.nst_randomizations}  \\
        --nst_distance         '${params.nst_distance}'      \\
        --nst_null_model       '${params.nst_null_model}'    \\
        ${aw_arg}                                            \\
        ${ses_arg}                                           \\
        ${rc_arg}                                            \\
        ${tax_arg}                                           \\
        ${excl_col}                                          \\
        ${excl_val}                                          \\
        ${params.group ? "--group '${params.group}'" : ""}
    """
}

// ---------------------------------------------------------------------------
// Process: QPE — betaNTI + Raup-Crick + EMS
// ---------------------------------------------------------------------------
process QPE {
    tag "${params.label}"

    publishDir "${params.output_dir}/QPE", mode: 'copy'

    input:
    path feature_table
    path meta_table
    path tree_file

    output:
    path "Coherence_*.csv",    emit: coherence
    path "Boundary_*.csv",     emit: boundary
    path "Turnover_*.csv",     emit: turnover
    path "Sitescores_*.csv",   emit: sitescores
    path "PairwiseRC_*.csv",   emit: pairwise_rc
    path "PairwisebNTI_*.csv", emit: pairwise_bnti
    path "RC_*.csv",           emit: rc
    path "bNTI_*.csv",         emit: bnti

    script:
    def tax_arg  = params.taxonomy_table ? "--taxonomy_table '${params.taxonomy_table}'" : ""
    def excl_col = params.exclude_column ? "--exclude_column '${params.exclude_column}'" : ""
    def excl_val = params.exclude_values ? "--exclude_values '${params.exclude_values}'" : ""

    """
    Rscript ${projectDir}/src/R/QPE.R \\
        --feature_table     '${feature_table}'             \\
        --meta_table        '${meta_table}'                \\
        --tree_file         '${tree_file}'                 \\
        --input_format      '${params.input_format}'       \\
        --output_dir        '.'                            \\
        --label             '${params.label}'              \\
        --scripts_dir       '${projectDir}'                \\
        --min_library_size  ${params.min_library_size}     \\
        --beta_reps         ${params.qpe_beta_reps}        \\
        --ems_sims          ${params.qpe_ems_sims}         \\
        ${tax_arg}                                         \\
        ${excl_col}                                        \\
        ${excl_val}                                        \\
        ${params.group ? "--group '${params.group}'" : ""}
    """
}

// ---------------------------------------------------------------------------
// Process: QPE_SUMMARY — assembly processes + EMS classification
// ---------------------------------------------------------------------------
process QPE_SUMMARY {
    tag "${params.label}"

    publishDir "${params.output_dir}/QPE", mode: 'copy'

    input:
    path pairwise_bnti
    path pairwise_rc
    path rc
    path coherence
    path boundary
    path turnover

    output:
    path "QPE_assembly_processes_*.csv", emit: assembly_processes
    path "EMS_*.csv",                   emit: ems

    script:
    def ordering_arg = params.qpe_ordering ? "--ordering '${params.qpe_ordering}'" : ""

    """
    Rscript ${projectDir}/src/R/QPE_summary.R \\
        --pairwise_bnti_csv '${pairwise_bnti}'        \\
        --pairwise_rc_csv   '${pairwise_rc}'          \\
        --rc_csv            '${rc}'                   \\
        --coherence_csv     '${coherence}'            \\
        --boundary_csv      '${boundary}'             \\
        --turnover_csv      '${turnover}'             \\
        --output_dir        '.'                       \\
        --label             '${params.label}'         \\
        ${ordering_arg}
    """
}

// ---------------------------------------------------------------------------
// Process: MERGE_PARQUET — merge all output CSVs into one Parquet file
// ---------------------------------------------------------------------------
process MERGE_PARQUET {
    tag "${params.label}"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path csvs

    output:
    path "environmental_filtering_${params.label}.parquet"

    script:
    """
    Rscript ${projectDir}/src/R/merge_parquet.R \\
        --label      '${params.label}' \\
        --output_dir '.'
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {

    // Feature table and metadata (always required)
    feat_ch = Channel.fromPath(params.feature_table, checkIfExists: true)
    meta_ch = Channel.fromPath(params.meta_table,    checkIfExists: true)

    // Phylogenetic tree (required by ENV_FILTERING and QPE)
    tree_ch = params.tree_file
        ? Channel.fromPath(params.tree_file, checkIfExists: true)
        : Channel.value(file("NO_FILE"))

    // Run analyses — NST runs independently of tree-based analyses
    ef_out  = ENV_FILTERING(feat_ch, meta_ch, tree_ch)
    nst_out = NST(feat_ch, meta_ch)
    qpe_out = QPE(feat_ch, meta_ch, tree_ch)

    // QPE_SUMMARY depends on QPE outputs
    sum_out = QPE_SUMMARY(
        qpe_out.pairwise_bnti,
        qpe_out.pairwise_rc,
        qpe_out.rc,
        qpe_out.coherence,
        qpe_out.boundary,
        qpe_out.turnover
    )

    // Optional: merge all CSVs into a single Parquet file
    if (params.merge_parquet) {
        all_csvs = ef_out.ef_wide
            .mix(ef_out.ef_long)
            .mix(ef_out.ef_pairwise)
            .mix(nst_out.nst_values)
            .mix(nst_out.nst_panova.ifEmpty(Channel.empty()))
            .mix(nst_out.nst_pairwise)
            .mix(qpe_out.coherence)
            .mix(qpe_out.boundary)
            .mix(qpe_out.turnover)
            .mix(qpe_out.sitescores)
            .mix(qpe_out.pairwise_rc)
            .mix(qpe_out.pairwise_bnti)
            .mix(qpe_out.rc)
            .mix(qpe_out.bnti)
            .mix(sum_out.assembly_processes)
            .mix(sum_out.ems)
            .collect()
        MERGE_PARQUET(all_csvs)
    }
}
