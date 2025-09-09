# Load necessary modules
module load StdEnv/2023
module load gcc/12.2.0 hmmer/3.3.2 blast+/2.14.0 augustus/3.5.0 prodigal/2.6.3


##BUSCO_Transcriptome completedness
busco -i trinity_out.Trinity.fasta  -l lineages -m transcriptome -o ./BUSCO