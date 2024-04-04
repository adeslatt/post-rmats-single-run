#!/bin/bash
#
#
#----------------------------------

allseqsaa=$(ls *.ORF_seqs_aa.fa)
seqsaa="_linear_aa.fa"
alllinearseq=$(ls *_linear_aa.fa)
blastdb="_blast.db"
blastresult="_blastp.tsv"

directory=$1
protein_query=$2

echo "directory is     " $directory
echo "protein_query is " $protein_query

for file in $alllinearseq; do

    name="${file%%_linear_aa.fa}"
    blastdbname=$name$blastdb
    blastoutname=$name$blastresult

    echo "name      is" $name
    echo "blastdbname is" $blastdbname
    echo "file        is" $file

    blastp -db $blastdbname -num_threads 8 -query $protein_query -gapopen 32767 -gapextend 32767 -outfmt "6 qseqid salltitles evalue length pident nident mismatch positive gapopen gaps bitscore qseq sseq" -max_target_seqs 1 -html -out $blastoutname
    
done


    
	  
