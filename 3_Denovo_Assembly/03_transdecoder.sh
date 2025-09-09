##Transdecoder
module load StdEnv/2020
module load gcc
module load transdecoder

TransDecoder.LongOrfs -t trinity_combined_cdhit.fasta
TransDecoder.Predict -t  trinity_combined_cdhit.fasta