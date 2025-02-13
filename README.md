# bulk-RNAseq-TREM2-Tgondii
Analysis files accompanying PLOS Pathogens paper entitled "TREM2 is a critical regulator of the immune response against acute Toxoplasma gondii infection."

This repository contains the Matlab code used to analyze and generate the images for the bulk-RNAseq analyses and figures from the paper. 

1. RNAseq_analysis_072624.m: this is a matlab code file that will run all the analyses and produce the figures from the paper. It was constructed in MATLAB R2024b version.
2. GO_bubble021325.m: a function for creating bubble plots using the list of differentially expressed genes

Major steps in analysis:
1. Pre-processing: Filtering out cells with low and high UMI, and high mitochondrial gene content.
2. Analysis of differentially expressed genes: Count data normalized to library size, then negative binomial exact test performed to determine the significance of observed expression differences between conditions.
3.  GO enrichment analysis: List of differentially expressed genes was inputted to g:GOSt platform on g:Profiler to undergo gene set enrichment analysis. Genes were mapped to data from the newest Ensembl databases to ensure that the matched terms are accurate.
