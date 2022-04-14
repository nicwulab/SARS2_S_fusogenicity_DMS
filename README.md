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

### Calculating experssion score and fusion score from DMS data
1. Merge overlapping paired-end reads using [PEAR](https://github.com/tseemann/PEAR)   
``pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]``
    - It can also be run using [./script/merge_reads.py](./script/merge_reads.py)
    - Output files should be placed in the fastq_merged/ folder and named as described in [./doc/filename_merged_fastq.tsv](./doc/filename_merged_fastq.tsv)

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
      - result/S2_HR1_DMS_count_aa.tsv

4. Compute expression score and fusion score   
``python3 script/S2_HR1_count2score.py``   
    - Input files:
      - result/S2_HR1_DMS_count_aa.tsv
    - Output file:
      - [./result/S2_HR1_DMS_scores.tsv](./result/S2_HR1_DMS_scores.tsv)
      - [./result/S2_HR1_DMS_scores_by_resi.tsv](./result/S2_HR1_DMS_scores_by_resi.tsv)
