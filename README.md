## Toxoplasma Gondii
Virulence of the parasite Toxoplasma Gondii depends on multiple factors that are packed to specific organelles. Although virulence factor expression is tightly regulated, the molecular mechanism controlling their regulation has been poorly understood. In this project we are taking the advantage of several computational approaches to study the regulators of such virulence, specifically, the Transcription Factors involved in lytic cycle of toxoplasma Gondii virulence traits.


### RNAseq Analysis Pipeline

***
This pipeline performs the following tasks:
* Performs quality control on FASTQ files using FASTQC
* Aligns the fatsq reads of each sample to the reference genome using Tophat
* Quantifies the expression of genes in each sample using cuffdiff from cufflinks package
* Find differentially expressed genes between various comparisons. 


