##Reference based assembly for Eukaryoytes (includes Splice variation)

##Step 1: Indexing 
module load star
STAR --runMode genomeGenerate --genomeDir Bd21-3_Index_GTF --genomeSAindexNbases 13 --genomeFastaFiles BdistachyonBd21_3_537_v1.0.fa --sjdbGTFfile Bd21-3.gtf
##Step 2: Alignment 
for i in *_R1_001_PE.fastq.gz; do
STAR  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts --sjdbGTFfile Bd21-3.gtf --genomeDir Bd21-3_Index_GTF --readFilesIn $i ${i%_R1_001_PE.fastq.gz}_R2_001_PE.fastq.gz --outFileNamePrefix ${i%_R1_001_PE.fastq.gz}/
done