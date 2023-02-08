## Studying the impact of HR1 mutations on SARS-CoV-2 spike expression and fusion using deep mutational scanning

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)

### Steps to analyze sequence conservation among variants of concern (VOC folder) and related sarbecoviruses (Sarbecovirus folder)
1. Download fasta files from genbank using Seq_downloader_1.0.py and GISAID (by hand).
2. Combine into 1 file using TxtProcesser.py.
3. Create a local database using makeblastdb (Local Blast, one for each alignment).
4. Use the reference sequence (SARS_CoV2_S_S2_ref_sequence.fasta) to run tblastn and generate BlastXML files (xx_BlastXML.xml).
5. Use XML_Extraction.py to extract information as csv files (xx_Blast_results.csv).
6. Use the extracted information to run MAFFT alignment (xx_Blast_results_alignment.txt).
7. Count residues at each position based on the reference sequence and calculate the alignment frequency (xx_sequence_conservation.csv).
8. Plot using xx_align_freq_vs_score.R.
