params.mode = 'full'

process PREP_ANCESTORS {
    publishDir "$params.resdir/$params.mode/ancestors", mode: 'copy'

    input:
    path input_data
    path region_info    // MAD4HATTER_Pf3D7_withPrimers.tsv
    path primers        // MAD4HATTER_Pf3D7_primers.tsv
    path code
    path sim_functions

    output:
    path "ancestors.rds", emit: ancestors

    script:
    """
    Rscript $code \
        --input           ${input_data} \
        --output_rds      ancestors.rds \
        --region_info     ${region_info} \
        --primers         ${primers} \
        --gen_pop_size    ${params.gen_pop_size} \
        --final_pop_size  ${params.final_pop_size} \
        --final_pop_gen2_prop ${params.final_pop_gen2_prop} \
        --anc_ancestor_size ${params.anc_ancestor_size} \
        --anc_pop_alpha       ${params.anc_pop_alpha} \
        --anc_coi_r           ${params.anc_coi_r} \
        --anc_coi_p           ${params.anc_coi_p} \
        --anc_k_s             ${params.anc_k_s} \
        --anc_max_k           ${params.anc_max_k} \
        --anc_max_coi         ${params.anc_max_coi}
    """
}

// ── Process 2: One simulation iteration (runs ×n_sims in parallel) ───────────
process RUN_SIM_ITER {
    publishDir "$params.resdir/$params.mode/iterations/iter_${iter_id}", mode: 'copy'

    input:
    path  ancestors           // ancestors.rds, staged into every work dir
    path  region_info       
    path  primers         
    path  code
    path  sim_functions
    val   iter_id             // e.g. 1..1000
    val   seed

    output:
    path  "iter_${iter_id}_dcifer.tsv.gz",  emit: dcifer_table
    path  "iter_${iter_id}_summary.tsv",    emit: summary

    script:
    """
    Rscript $code \
        --ancestors_rds     ${ancestors} \
        --region_info       ${region_info} \
        --primers           ${primers} \
        --iter_id           ${iter_id} \
        --seed              ${seed} \
        --sim_pop_size      ${params.sim_pop_size} \
        --num_samples       ${params.num_samples_per_sim} \
        --pop_alpha         ${params.pop_alpha} \
        --coi_r             ${params.coi_r} \
        --coi_p             ${params.coi_p} \
        --k_s               ${params.k_s} \
        --max_k             ${params.max_k} \
        --max_coi           ${params.max_coi} \
        --dcifer_output     iter_${iter_id}_dcifer.tsv.gz \
        --summary_output    iter_${iter_id}_summary.tsv
    """
}

// ── Process 3: Combine all summary rows into one table ────────────────────────
process COMBINE_SUMMARIES {
    publishDir "$params.resdir/$params.mode", mode: 'copy'

    input:
    path summaries    // all *_summary.tsv files collected into the work dir

    output:
    path "all_summaries.tsv"

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)

    files <- list.files(".", pattern = "_summary\\\\.tsv\$", full.names = TRUE)
    combined <- purrr::map_dfr(files, read_tsv, col_types = cols())
    write_tsv(combined, "all_summaries.tsv")
    """
}


workflow {
    // One shared seed for ancestor prep; per-sim seeds derived from iter_id
    // so results are fully reproducible
    ancestor_seed = Channel.value(12345)

    // Step 1: build ancestors once from the shared input
    PREP_ANCESTORS(
        file(params.input),
        file('data/MAD4HATTER_Pf3D7_withPrimers.tsv'),
        file('data/MAD4HATTER_Pf3D7_primers.tsv'),
        file('bin/scripts/prep_ancestors.R'),
        file('bin/scripts/simulation_functions.R')
    )

    // Step 2: create a channel of 1..n_sims iteration IDs
    // In test mode, just run the first 3
    n = (params.mode == 'test') ? 3 : params.n_sims

    iter_ids = Channel.from(1..n)

    // Combine the single ancestors file with every iter_id
    // PREP_ANCESTORS.out.ancestors is a value channel after .first(),
    // so it broadcasts to all iter_ids without being consumed
    RUN_SIM_ITER(
        PREP_ANCESTORS.out.ancestors,
        file("$projectDir/data/MAD4HATTER_Pf3D7_withPrimers.tsv"),  // <-- add
        file("$projectDir/data/MAD4HATTER_Pf3D7_primers.tsv"),
        file('bin/scripts/run_sim_iter.R'),
        file('bin/scripts/simulation_functions.R'),
        iter_ids,
        iter_ids.map { id -> id + 99999 }   // unique but reproducible seed per iter
    )

    // Step 3: wait for all summaries then combine
    COMBINE_SUMMARIES(
        RUN_SIM_ITER.out.summary.collect()
    )
}
