module load subread
featureCounts -a Bd21-3.gtf -p --countReadPairs -s 2 -t exon -g gene_id -o ./Feature.count_Stranded.txt ./*.bam