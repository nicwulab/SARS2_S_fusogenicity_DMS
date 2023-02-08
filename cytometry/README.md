## Studying the impact of HR1 mutations on SARS-CoV-2 spike expression and fusion using deep mutational scanning: Flow cytometry validation

### Dependencies
* [R](https://www.r-project.org/) (version 4.1)

### Input files
* [./expr.csv](./expr.csv): Median fluorescence intensity (MFI) values for assaying surface expression of membrane-bound spike mutants
* [./fus3.csv](./fus3.csv): Median fluorescence intensity (MFI) values for assaying fusion of membrane-bound spike mutants with hACE2-expressing cells at 3 hours post-mixing

### Generating output
* Run the R script to generate plot of MFI values relative to that of WT for membrane-bound spike mutations

### Notes
* 'Pre', 'post', and 'combo' refer to candidate prefusion-stabilizing mutations, fusion-competent mutations, and combinations of candidate prefusion-stabilizing mutations, respectively.
