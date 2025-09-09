module load hisat2
##Step 1: Indexing
hisat2-build --ss Bd21-3_splicesites.tsv --exon Bd21-3_exons.tsv BdistachyonBd21_3_537_v1.0.fa Bd21-3_Index_HiSAT2
##Step 2: Alignment 
hisat2 -x ./Bd21-3_Index_HiSAT2 --rna-strandness RF -1 S0_L00_R1_001_PE.fastq.gz -2 S0_L00_R2_001_PE.fastq.gz -S Sample001.bam
