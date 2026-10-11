process SALMON_INDEX {
    tag "Building Decoy-Aware Index"
    
    publishDir "${params.outdir}/reference/salmon_index", mode: 'copy'

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "true_salmon_index", emit: index_dir

    script:
    """
    # 1. Isolate structural chromosomes to serve as decoys
    grep "^>" ${genome_fasta} | cut -d ' ' -f 1 > decoys.txt
    sed -i -e 's/>//' decoys.txt
    
    # 2. Combine files into a clean master sequence file
    cat ${transcript_fasta} ${genome_fasta} > gentrome.fa
    
    # 3. Compile the pristine index
    salmon index -t gentrome.fa -d decoys.txt -p ${task.cpus} -i true_salmon_index
    """
}