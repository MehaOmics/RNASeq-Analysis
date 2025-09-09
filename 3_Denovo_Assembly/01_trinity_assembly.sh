##Denovo based
module load StdEnv/2023
module load gcc/12.3
module load salmon/1.10.2
module load openmpi/4.0.3
module load trinity/2.14.0
module load samtools
module load jellyfish

Trinity --seqType fq  --max_memory 200G --SS_lib_type RF --left all_reads_R1.fastq.gz   --right all_reads_R2.fastq.gz --CPU 40 --no_bowtie --output ./Trinity