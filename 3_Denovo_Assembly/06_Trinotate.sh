# Load Trinotate module
module load StdEnv/2020 gcc/9.3.0 trinotate/4.0.0


# Step 1: Create a Trinotate SQLite database
Trinotate --db Trinotate.sqlite --create --trinotate_data_dir trinotate_database



# Step 2: Initialize the database with transcripts and predicted proteins
Trinotate --db Trinotate.sqlite --init --gene_trans_map trinity_out.Trinity.fasta.gene_trans_map --transcript_fasta trinity_out.Trinity.fasta --transdecoder_pep trinity_out.Trinity.fasta.transdecoder.pep 

#To check if  database has loaded correctly:
sqlite3 Trinotate.sqlite ".tables"

## Trinotate requires results of blastx, blastp, signalp, pfam, tmhmm, eggNOG and cmscan  

##Run blastx and blastp, your data as query and use the uniprot database downloaded above
## For BLAST: Make database for uniprot

#makeblastdb -in uniprot_sprot.fasta -dbtype prot -out my_blast_db

##Run blastp for homology based seqeunce search 
blastp -query trinity_out.Trinity.fasta.transdecoder.pep -db my_blast_db -num_threads 10 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

##Run blastx for homology based seqeunce search 
blastx -query trinity_out.Trinity.fasta -db my_blast_db -num_threads 10 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

##RunHMMScan for homology based seqeunce search 
hmmscan --domtblout pfam.out Pfam-A.hmm trinity_out.Trinity.fasta.transdecoder.pep

##Run signalp for Computational prediction of sequence features
signalp6 --fastafile trinity_out.Trinity.fasta.transdecoder.pep --organism euk --output_dir ./ --format txt --mode fast


##Run TMHMM to predict transmembrane regions
tmhmm --short < trinity_out.Trinity.fasta.transdecoder.pep > tmhmm.out

#Download eggNOG database
python download_eggnog_data.py --data_dir ./eggnog-mapper

#RUNNING EGGnog
emapper.py -i ../trinity_out.Trinity.fasta.transdecoder.pep -o ./output/EggNOG_Pep --cpu 20  --data_dir ./

##Run cmscan
cmscan --cpu 10 --tblout rfam.out --fmt 2 --nohmmonly Rfam.cm trinity_out.Trinity.fasta > rfam.cmscan.out

####Now load all the results in Trinotate database (Trinotate.sqlite)
# Load BLAST results (SwissProt)
Trinotate --db Trinotate.sqlite --LOAD_swissprot_blastx blastx.outfmt6
Trinotate --db Trinotate.sqlite --LOAD_swissprot_blastp blastp.outfmt6

# Load Pfam domain hits (from hmmscan)
Trinotate --db Trinotate.sqlite --LOAD_pfam pfam.out

# Load SignalP predictions
Trinotate --db Trinotate.sqlite --LOAD_signalp6 signalp.out

# Load TMHMM predictions
Trinotate --db Trinotate.sqlite --LOAD_tmhmmv2 tmhmm.out

# Load EggNOG-mapper annotations
Trinotate --db Trinotate.sqlite --LOAD_EggnogMapper emapper.annotations

# load Infernal results for structural RNA 
Trinotate --db Trinotate.sqlite --LOAD_infernal cmscan.out

# Generate final annotation report
Trinotate --db Trinotate.sqlite --report > trinotate_annotation_report.xls

#Summary stats

perl /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/trinotate/4.0.0/util/report_summary/trinotate_report_summary.pl trinotate_annotation_report.xls Trinotate_report_stats



