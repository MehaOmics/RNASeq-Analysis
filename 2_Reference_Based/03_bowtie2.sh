##Reference based assembly for Bacteria (using bowtie, no splice variation)
module load bowtie2
module load samtools

##Step 1: Indexing
bowtie2-build -f B26_genome.fasta bow2_index_B26

##Step 2: Alignment
bowtie2  -x bow2_index_B26 -1 ../S0_L001_R1_001_PE.fastq.gz -2 ../S0_L001_R2_001_PE.fastq.gz  --very-sensitive | samtools view -bS - > Sample.bam