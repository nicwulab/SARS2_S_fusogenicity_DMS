## Studying the impact of HR1 mutations on SARS-CoV-2 spike expression and fusion using deep mutational scanning

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR)
* [ProDy](http://prody.csb.pitt.edu/)

### Input files
* [./Fasta/HR1_ref_seq.fa](./Fasta/HR1_ref_seq.fa): HR1 WT amino acid sequence
* [./Fasta/sars-cov-2-spike-hr1-linker.fa](./Fasta/sars-cov-2-spike-hr1-linker.fa): HR1 nucleotide sequence with 21 nt upstream (5' flank)
* Raw read files in fastq format from NIH SRA database [BioProject PRJNA826665](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA826665)

### Primer design for DMS library construction
1. Generating foward (NNK + internal barcode) and reverse primers (constant)   
``python3 script/lib_primer_design.py``
    - Input file:
      - [./Fasta/sars-cov-2-spike-hr1-linker.fa](./Fasta/sars-cov-2-spike-hr1-linker.fa)
    - Output files:
      - [./Fasta/SARS-2-Spike-HR1_lib_Fprimer_bc.fa](./Fasta/SARS-2-Spike-HR1_lib_Fprimer_bc.fa)
      - [./Fasta/SARS-2-Spike-HR1_lib_Rprimer_bc.fa](./Fasta/SARS-2-Spike-HR1_lib_Rprimer_bc.fa)

2. Generating barcode file   
``python3 script/check_barcode.py``
    - Input files:
      - [./Fasta/SARS-2-Spike-HR1_lib_Fprimer_bc.fa](./Fasta/SARS-2-Spike-HR1_lib_Fprimer_bc.fa)
      - [./Fasta/sars-cov-2-spike-hr1-linker.fa](./Fasta/sars-cov-2-spike-hr1-linker.fa)
    - Output file:
      - [./data/barcodes.tsv](./data/barcodes.tsv)

### I. Calculating experssion score and fusion score from DMS data
1. Merge overlapping paired-end reads using [PEAR](https://github.com/tseemann/PEAR)   
``pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]``
    - It can also be run using [./script/merge_reads.py](./script/merge_reads.py)
    - Output files should be placed in the fastq_merged/ folder and named as described in [./data/filename_merged_fastq.tsv](./data/filename_merged_fastq.tsv)

2. Counting variants based on nucleotide sequences   
``python3 script/S2_HR1_fastq2count.py``   
    - Input files:
      - Merged read files in fastq_merged/ folder
    - Output file:
      - result/S2_HR1_DMS_count_nuc.tsv

3. Convert nucleotide sequences to amino acid mutations   
``python3 script/S2_HR1_count_nuc2aa.py``   
    - Input file:
      - [./data/barcodes.tsv](./data/barcodes.tsv)
      - [./Fasta/HR1_ref_seq.fa](./Fasta/HR1_ref_seq.fa)
      - result/S2_HR1_DMS_count_nuc.tsv
    - Output file:
      - [./result/S2_HR1_DMS_count_aa.tsv](./result/S2_HR1_DMS_count_aa.tsv)

4. Compute expression score and fusion score   
``python3 script/S2_HR1_count2score.py``   
    - Input files:
      - [./result/S2_HR1_DMS_count_aa.tsv](./result/S2_HR1_DMS_count_aa.tsv)
    - Output file:
      - [./result/S2_HR1_DMS_scores.tsv](./result/S2_HR1_DMS_scores.tsv)
      - [./result/S2_HR1_DMS_scores_by_resi.tsv](./result/S2_HR1_DMS_scores_by_resi.tsv)

### II. Checking data quality 
1. Plot correlation between replicates and compare silent/missense/nonsense   
``Rscript script/plot_QC.R``
    - Input file:
      - [./result/S2_HR1_DMS_scores.tsv](./result/S2_HR1_DMS_scores.tsv)
    - Output files:
      - [./graph/QC_replicate_exp12.png](./graph/QC_replicate_exp12.png)
      - [./graph/QC_replicate_exp13.png](./graph/QC_replicate_exp13.png)
      - [./graph/QC_replicate_exp23.png](./graph/QC_replicate_exp23.png)
      - [./graph/QC_replicate_fus.png](./graph/QC_replicate_fus.png)
      - [./graph/Exp_by_class.png](./graph/Exp_by_class.png)
      - [./graph/Fus_by_class.png](./graph/Fus_by_class.png)

2. Plot heatmap for the expression scores of individual mutations   
``Rscript script/plot_score_heatmap``
    - Input file:
      - [./result/S2_HR1_DMS_scores.tsv](./result/S2_HR1_DMS_scores.tsv)
    - Ouput file:
      - [./graph/S2_HR1_fus_heatmap.png](./graph/S2_HR1_fus_heatmap.png)

### III. Identify prefusion-stabilizing mutations
1. Plot the relationship between fusion score and expression score
``script/plot_exp_vs_fus.R``
    - Input file:
      - [./result/S2_HR1_DMS_scores.tsv](./result/S2_HR1_DMS_scores.tsv)
    - Output file: 
      - [./graph/Exp_vs_fus.png](./graph/Exp_vs_fus.png)
      - [./graph/Exp_vs_fus_reidual.png](./graph/Exp_vs_fus_reidual.png)
