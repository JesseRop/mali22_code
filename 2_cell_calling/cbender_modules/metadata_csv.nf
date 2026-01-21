process METADATA_CSV {
    tag { h5.getName() }
    cpus 1
    memory '20 GB'
    time '1.h'

    publishDir { o_dir }, 
        mode: 'copy', 
        pattern: "*.csv", 
        overwrite: true

    input:
        tuple val(msc_id), val(cbender), path(h5), val(o_dir)

    output:
        tuple path("${params.obs_fname}"), path("${params.var_fname}"), emit: metadata_csvs

    script:
    """
    set -euo pipefail

    source ${params.cellbender_env}/bin/activate
    
    python ${params.metadata_script} \\
        --input "${h5}" \\
        --outdir . \\
        --obs-fname "${params.obs_fname}" \\
        --var-fname "${params.var_fname}"
    """
}