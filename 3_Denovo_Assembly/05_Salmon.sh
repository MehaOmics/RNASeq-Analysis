module load salmon

##Step 1: Salmon index
salmon index -t trinity_out.Trinity.fasta -i Salmon_Index -k 31

##Step 2: Quantification
salmon quant -i Salmon_Index -l A -1 ../PE_R1.fastq.gz -2 ../PE_R2.fastq.gz  --validateMappings -o ./Salmon_quant/S1

