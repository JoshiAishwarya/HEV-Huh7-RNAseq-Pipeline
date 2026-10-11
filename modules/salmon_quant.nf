process SALMON_QUANT {
    tag "$sample_id"
    publishDir "${params.outdir}/salmon_quant", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index_dir

    output:
    path "${sample_id}", emit: out_dir
    path "${sample_id}/quant.sf", emit: quant

    script:
    """
    salmon quant \\
        -i ${index_dir} \\
        -l A \\
        -r ${reads} \\
        -p ${task.cpus} \\
        --validateMappings \\
        -o ${sample_id}
    """
}