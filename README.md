# HEV Exosome-Pathway RNA-seq Quantification & Differential Expression Pipeline

A production-grade Nextflow pipeline for end-to-end gene-level expression quantification and downstream statistical analysis from Ion Torrent single-end RNA-seq data (HEV replicon-transfected Huh7 cells, PRJNA381374). 

The pipeline handles raw read quality preprocessing, implements a strict **decoy-aware indexing framework** to filter genomic false-positives, and is paired with a multi-factor **DESeq2 interaction model** in R/RStudio for high-accuracy differential expression testing and functional pathway annotation.

## Pipeline Architecture





```
Input FASTQ Files (Raw_data/*.fastq.gz)
│
├──> FASTQC_RAW (Initial Read Quality Control)
│
└──> FASTP (Adapter Trimming & Low-Phred Filtering)
│
├──> Weightless Zero-Byte Symlink Staging (stageInMode = 'symlink')
│
└──> SALMON_QUANT (Transcript Abundance Profiling via Decoy Index)
│
├──> MULTIQC (Aggregated Quality Assessment Dashboard)
└──> Results/salmon_quant//quant.sf
│
└──> [Downstream R/RStudio Statistical Interface Engine]
├──> tximport (Decimal-stripped Transcript Aggregation)
├──> DESeq2 (Multi-Factor Interaction Model Matrix)
├──> org.Hs.eg.db (Official Gene Symbol Annotation)
└──> clusterProfiler (Genome-Wide Functional Enrichment)

```


## Reference Optimization & Decoy-Aware Indexing

To prevent false-positive expression spikes, this pipeline implements a modern **decoy-aware reference framework** instead of a naked transcriptome alignment loop. 

Naked transcript files (`cdna.fa`) create a mathematical blind spot where fragments originating from unannotated genomic regions or repetitive DNA elements are erroneously assigned to transcripts, skewing log2 fold changes. This pipeline dynamically concatenates **Ensembl's GRCh38 cDNA transcripts** with full **chromosomal genomic reference sequences (`*_genomic.fna`)** to construct a true biological decoy index. Genomic background noise reads are successfully absorbed by the decoy chromosomes, leaving only authentic, unconfounded expression metrics in your final dataset.

## Core Feature Enhancements
* **Zero-Disk Overhead Staging:** Configured with `stageInMode = 'symlink'` to process hundreds of gigabytes of sequence data through zero-byte virtual shortcuts without duplication.
* **Automated Workspace Garbage Cleanup:** Enforces an internal `cleanup = true` rule to automatically purge intermediate execution tracking logs upon success, protecting local disk capacity.

## Tools Ecosystem
* **FastQC:** Automated quality metrics evaluation of raw sequences.
* **fastp:** High-performance adapter and quality trimming (Phred ≥ 20, minimum length filtering ≥ 36bp).
* **Salmon:** Decoy-aware pseudo-alignment k-mer indexing and quantification. Generates length-corrected abundance profiles (`quant.sf`) without heavy BAM disk utilization.
* **MultiQC:** Consolidates raw tool execution logs into a centralized dashboard report.
* **Nextflow:** Workflow parallelization orchestration, algorithmic caching, and strict check-pointed resume tracking.
* **DESeq2:** Adaptive negative binomial modeling built for small-sample, multi-factor variance estimation.
* **clusterProfiler:** Programmatic over-representation profiling engine mapped to the Gene Ontology database.


## Installation & Setup

```bash
# Clone the repository
git clone https://github.com
cd Nextflow-dev

# Initialize the localized bioinformatic environment
conda env create -f environment.yml
conda activate bnf
```

## Core Execution Commands

Execute the pipeline across local storage bounds using default parameter pathways:
```bash
nextflow run main.nf
```

Resume an execution stream seamlessly from its last healthy checkpoint:
```bash
nextflow run main.nf -resume
```

## Downstream Statistical Analysis (DESeq2 in RStudio)

Once Nextflow successfully generates your clean quantification files (`Results/salmon_quant/`), downstream statistical parsing is handled natively in R/RStudio to factor out timeline batch effects and map clear biological trajectories.

### The Interaction Design Matrix
Your sample metrics span two distinct timeline checkpoints (24h acute phase vs 72h late phase). Subsetting these time blocks into isolated matrices breaks statistical power and induces a zero-variance artifact because the 72h phase contains only a single baseline control replicate (`SRR5407648`). 

To resolve this, the analysis uses all 7 samples concurrently inside a unified **Interaction Model Design**:
\[\text{design} = \sim \text{time} + \text{condition} + \text{time:condition}\]

This architecture pools global variance across all controls simultaneously, stabilizes statistical confidence, suppresses false-positive zero-inflation spikes, and lowers log2 fold changes to realistic biological parameters (\(\pm1\) to \(\pm4\)).

### Automated Annotation & Subsetting
The R interface strips transcript decimals (`.6`, `.3`) on the fly for a 1:1 match with the **`org.Hs.eg.db`** annotation database, converting obscure `NM_` accession headers into **Official Gene Symbols** (e.g., *STAT1*, *IFI44*). It extracts extreme significant responders (\(padj < 0.05\)) across **6 clean summary spreadsheets**:
* `time24h_top20_upregulated.csv` & `time24h_top20_downregulated.csv` (Early immune defenses)
* `time72h_top20_upregulated.csv` & `time72h_top20_downregulated.csv` (Late viral hijack patterns)
* `overall_top20_upregulated.csv` & `overall_top20_downregulated.csv` (Core longitudinal variants)

### Functional Over-Representation Enrichment
To extract biological meaning from the full master tables, the pipeline passes all genes meeting \(padj < 0.05\) through **`clusterProfiler::enrichGO`** to evaluate **Molecular Functions (MF)**. This yields independent, genome-wide pathway charts mapping out the kinetic switch of the virus:
* **24h Interaction Profile:** Showcases significant activation of `single-stranded RNA binding` and `helicase activity`, capturing active viral genome uncoiling.
* **72h Interaction Profile:** Showcases a shift to `protein serine/threonine kinase activity`, capturing downstream metabolic stress and cell signaling events.

## Pipeline Deliverables Dashboard

All final production deliverables are neatly exported to your localized target output path:

### Sequence Processing Directory
* `Results/salmon_quant/<sample_id>/quant.sf` — Length-corrected transcript abundances.
* `Results/multiqc/multiqc_report.html` — Centralized raw read processing quality dashboard.

### Statistical & Pathway Analysis Directory
* `deseq2_24h_results.csv` & `deseq2_72h_results.csv` — Full annotated master files.
* `time24h_top20_upregulated.csv` & `time72h_top20_downregulated.csv` — Subset summaries.
* `pca_plot.png` — Strong spatial sample clustering by condition and time boundaries.
* `expression_heatmap.png` — Hierarchical clustering of the top 50 highly variable genes labeled with official symbols.
* `volcano_24h.png` & `volcano_72h.png` — Symmetrical volcano distributions with matching red highlights.