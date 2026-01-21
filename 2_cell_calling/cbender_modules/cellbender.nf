process CELLBENDER {
    tag { "${msc_id} - LR${l_rate}_E${epoch} ${cbender}" }
    cpus { params.ncores }
    memory { params.mem + ' MB' }
    time '5.h'
    queue 'gpu-normal'
    // errorStrategy 'ignore'
    clusterOptions "-gpu 'num=1:j_exclusive=yes'"

    publishDir { "${o_dir}" }, 
        mode: 'copy', 
        pattern: "cb*", 
        overwrite: true

    input:
        tuple val(msc_id), val(irods_id), val(cbender), path(raw_mtx), 
              val(o_dir), val(o_file), val(l_rate), val(epoch), 
              path(exp_cells_file), path(exp_edrops_file)

    output:
        tuple val(msc_id), val(cbender), path("${o_file}"), val(o_dir), path("cb*"), emit: cb_outputs

    script:
    """
    mkdir -p ${o_dir}
    
    exp_cells=\$(cat ${exp_cells_file})
    exp_edrops=\$(cat ${exp_edrops_file})

    ${params.cellbender_script} ${raw_mtx} ${o_file} ${l_rate} ${epoch} \$exp_cells \$exp_edrops 
    """
} 