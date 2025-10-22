# 2025-dibrotation-snakemake
Snakemake workflow to perform metagenome analysis on single sample

## Metagenomics Rotation Project
This workflow is the result from a project to become comfortable with bioinformatics tools. The [project](https://dib-lab.github.io/dib_rotation/) comes from the DIB lab at UC Davis and utilizes data reported in a [2016 paper by Hu et al](https://journals.asm.org/doi/10.1128/mbio.01669-15).  
  
### Summary
In short, this workflow takes reads from a single sample and performs various rules:
* Quality assessment with fastqc
* Quality trimming with fastp
* k-mer trimming with khmer
* Separate comparisons of trimmed and raw data with sourmash `sketch` and `gather`
* Bin completion with spacegraphcats
  * Important: choice of bin is from a priori analysis
  
### Resulting File Structure
