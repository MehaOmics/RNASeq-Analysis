## RNA-seq analysis
- This repository is a compilation of scripts used during my PhD project to study transcriptional changes during plant interactions with plant growth promoting bacteria
- It includes scripts for **RNA-seq analysis** for both **reference-based** and **de novo** transcriptome analysis  
- The repository contains scripts for alignment, assembly, annotation, quantification and visualization using R


## Folder Structure
1. Quality_control: Scripts for Quality check using Fastqc and trimming of adapters (Trimmomatic)
2. Reference_Based: Scripts for reference-based workflows (STAR, HISAT2, Bowtie2 and FeatureCounts)
3. Denovo_Assembly: Scripts for de novo assembly, quantification and annotation (Trinity, CD-HIT, TransDecoder, BUSCO, Salmon and Trinnotate)
4. R_Analysis: R scripts for Diffrential Gene Expression analysis (DeSeq2, EdgeR) and R scripts for generating PCA plots, heatmaps, and volcano plots

## Requirements
- Linux/HPC environment (tested on Compute Canada)
- Tools: FastQC, Trimmomatic, STAR, HISAT2, Bowtie2, Trinity, CD-HIT, TransDecoder, BUSCO, Salmon
- R with packages: ggplot2, ComplexHeatmap, pheatmap

## Author
- Meha Sharma
