include { FASTQC as FASTQC_RAW }     from '../modules/fastqc.nf'
include { FASTP }                    from '../modules/fastp.nf'
include { SALMON_INDEX }             from '../modules/salmon_index.nf'
include { SALMON_QUANT }             from '../modules/salmon_quant.nf'
include { MULTIQC }                  from '../modules/multiqc.nf'

workflow RNASEQ {
    take:
    read_ch
    genome_fasta
    transcript_fasta

    main:
    // 1. Quality Control on RAW read profiles
    FASTQC_RAW(read_ch)
    
    // 2. Smart data trimming via fastp
    FASTP(read_ch)

    // 3. Compile decoy-aware Salmon Reference Index 
    idx_dir = SALMON_INDEX(genome_fasta, transcript_fasta)
    
    // 4. Calculate high-accuracy target expression values
    quant_out = SALMON_QUANT(FASTP.out.trimmed_reads, idx_dir)

    // 5. Aggregate operational runs across MultiQC Dashboard
    mqc_inputs = FASTQC_RAW.out
        .mix(FASTP.out.json, quant_out.out_dir)
        .collect()
    MULTIQC(mqc_inputs)
}