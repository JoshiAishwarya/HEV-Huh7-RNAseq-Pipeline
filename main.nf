#!/usr/bin/env nextflow

include { RNASEQ } from './workflows/rnaseq.nf'

workflow {
    log.info """\
        HEV Huh7 RNA-seq Production Pipeline (Core Sequence Processing)
        ==================================================================
        Reads directory  : ${params.reads}
        Reference Genome : ${params.genome_dir}
        Local Output Path: ${params.outdir}
        """.stripIndent()

    // Initialize fastq sample data streams
    read_ch = Channel
        .fromPath(params.reads, checkIfExists: true)
        .map { f -> tuple(f.simpleName, f) }

    // Locate genomics files dynamically inside the target folder
    genome_fasta_ch     = Channel.fromPath("${params.genome_dir}/*_genomic.fna", checkIfExists: true).collect()
    transcript_fasta_ch = Channel.fromPath("${params.genome_dir}/rna.fna", checkIfExists: true).collect()

    // Kick off core processing with exactly 3 matching arguments
    RNASEQ(read_ch, genome_fasta_ch, transcript_fasta_ch)
}