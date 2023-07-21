#!/bin/bash
#---------------------------------------------------------------------
#
# name:  extract_SE_coordinates.sh
# purpose: take the output from rMATS with the coordinates of each
#          splice defined
#          1. extract the SE coordinates IJC and SJC numbers
#             GeneID
#             geneSymbol
#             chr
#             strand
#             exonStart_0base
#             exonEnd
#             upstreamEE
#             downstreamES
#             IJC_SAMPLE_1
#             SJC_SAMPLE_2
# input:   the directory
# output:  the coordinates and gene identifiers, strand for SE
#---------------------------------------------------------------------

cd $1
allSE="*.SE.MATS.JC.txt"
SEend=".SE.txt"
PWD=$(pwd)
echo "Current Working Directory is = " $PWD
echo "allSE                        = " $allSE


# loop through 
for file in $allSE; do
    name="${file%-*.SE.MATS.JC.txt}"
    name_se=$name$SEend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_se             = " $name_se
    
    cut -f 1-14 $file > $name_se
done

