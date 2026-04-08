params.mode = 'full'

process DCIFER_WRAPPER_PER_POP_AF {
    publishDir "$params.resdir/$params.mode/rel/$sim_id", mode: 'copy'

    def btwn_host_rel_output = 'btwn_host_rel_output.tsv'

    input:
    path code
    tuple val(sim_id), path(allele_table)
    // Allows passing of arbitrary additional arguments without 
    // creating a different process
    val extra_args
    val seed

    output:
    tuple val(sim_id), path("${btwn_host_rel_output}.gz"), emit: btwn_host_rel

    script:
    """
    Rscript $code \
        --allele_table ${allele_table} \
        $extra_args \
        --threads ${task.cpus} \
        --seed $seed \
        --btwn_host_rel_output ${btwn_host_rel_output}
    gzip ${btwn_host_rel_output}
    """
}

workflow {

    // Workflow inputs
    seed = Channel.value(25079987)
    // This creates a channel for all of your allele tables. A channel 
    // is a queue that contains all of the inputs you want to process 
    // (in this case, allele tables to estimate relatedness from).
    allele_tables = channel
        .fromPath("$params.simdir/*/allele_table.tsv")
    // If we are running in test mode, just use the first allele table
    if (params.mode == 'test') {
        allele_tables = allele_tables.first()
    }

    // This block of code converts the channel of files into a channel 
    // of tuples, where each contains the simulation ID and the file 
    // path. This is useful later for keeping results organized.
    allele_tables = allele_tables
        .map { file -> 
            // Get name of parent directory
            def sim_id = file.parent.name
                // Strip off "sim_"
                .replaceAll(/sim_$/, '')
            return [sim_id, file]
        }

    // Run Dcifer on all allele tables
    DCIFER_WRAPPER_PER_POP_AF(
        // The code is included as an input so that Nextflow will 
        // automatically rerun everything if the code changes. If you 
        // don't want that, just delete this from here and the process.
        file('bin/scripts/dcifer_wrapper_per_pop_af.R'), 
        allele_tables, 
        '', 
        seed
    )

}
