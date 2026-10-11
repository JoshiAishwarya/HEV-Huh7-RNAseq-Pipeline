process FASTP {
    tag "$sample_id"
    
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*.trimmed.fastq.gz"
    publishDir "${params.outdir}/fastp_reports", mode: 'copy', pattern: "*.{json,html}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed_reads
    path("${sample_id}.fastp.json"), emit: json
    path("${sample_id}.fastp.html")

    script:
    """
    fastp \\
        -i ${reads} \\
        -o ${sample_id}.trimmed.fastq.gz \\
        --qualified_quality_phred ${params.qual_cutoff} \\
        --length_required ${params.min_len} \\
        --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html
    """
}