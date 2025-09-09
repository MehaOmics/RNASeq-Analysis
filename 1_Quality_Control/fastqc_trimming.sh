##Step 1: Quality check for raw fasta files
module load fastqc
fastqc *.gz --outdir=./fqc_reports/

##Step 2: Trimming of adapters (For loop for all the fastq files in the folder)

for R1 in *R1*
do
   R2=${R1//R1_001.fastq.gz/R2_001.fastq.gz}
   R1PE=${R1//.fastq.gz/_PE.fastq.gz}
   R1SE=${R1//.fastq.gz/_SE.fastq.gz}	
   R2PE=${R2//.fastq.gz/_PE.fastq.gz}
   R2SE=${R2//.fastq.gz/_SE.fastq.gz}	
   java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 $R1 $R2 $R1PE $R1SE $R2PE $R2SE ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50
done