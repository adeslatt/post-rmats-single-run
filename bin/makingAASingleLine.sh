#!/bin/bash
#
#
#----------------------------------

allseqsaa=$(ls *.ORF_seqs_aa.fa)
seqsaa="_linear_aa.fa"


for file in $allseqsaa; do

    name="${file%%.ORF_seqs_aa.fa}"
    outname=$name$seqsaa

    echo "name is   " $name
    echo "outname is" $outname
    echo "file    is" $file
    
    awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' "$file" > $outname

done


    
	  
