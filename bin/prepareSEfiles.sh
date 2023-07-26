#!/bin/bash
#---------------------------------------------------------------------
#
# name:  extract_SE_coordinates.sh
#
# philosophy - always an elements of style approach
#
#  input - process - output
#
#
# input:   the directory
#
# process:
#
#   step 1 - extract SE coordinages with IJC and SJC
#
#          take the output from rMATS with the coordinates of each
#          splice defined
#          extract the SE coordinates IJC and SJC numbers
#          col 1 - ID -- DO NOT NEED/USE
#          col 2 - GeneID
#          col 3 - geneSymbol
#          col 4 - chr
#          col 5 - strand
#          col 6 - exonStart_0base
#          col 7 - exonEnd
#          col 8 - upstreamES
#          col 9 - upstreamEE
#          col 10 - downstreamES
#          col 11 - downstreamEE
#          col 12 - ID - DO NOT NEED/USE
#          col 13 - IJC_SAMPLE_1
#          col 14 - SJC_SAMPLE_2
#
#          remove the header
#
#   step 2 - cat the SE files with coordinates and counts together (create a union)
#         
#   step 3 - unique sort master union file (all.SE.txt)
#
#   step 4 - cut the last two columns from sorted unique union file
#
#   step 5 - add an ID
#
#   step 6 - sort the individual files
#
#   step 7 - normalize the individual files
#
# output:  the coordinates and gene identifiers, strand for SE
#---------------------------------------------------------------------

cd $1
allSE="*.SE.MATS.JC.txt"
SEend=".SE.txt"
tmp="tmp.txt"
PWD=$(pwd)
echo "Current Working Directory is = " $PWD
echo "allSE                        = " $allSE

#
#   step 1 - extract coordinates from the SE file
#
for file in $allSE; do
    name="${file%-*.SE.MATS.JC.txt}"
    name_se=$name$SEend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_se             = " $name_se
    
    cut -f 2-11,13-14 $file > $name_se

    # remove the header
    tail -n +2 $name_se > $tmp && mv $tmp $name_se
done

#
#   step 2 - create the master union file by cat'ing these files together
#
allSEend="*.SE.txt"
allSE="SE.all.txt"

cat $allSEend > $allSE

#
#   step 3 - unique sort master union file (all.SE.txt)
#
# sort unique
#   our column structure has changed
#          col 1 - GeneID
#          col 2 - geneSymbol
#          col 3 - chr
#          col 4 - strand
#          col 5 - exonStart_0base
#          col 6 - exonEnd
#          col 7 - upstreamES
#          col 8 - upstreamEE
#          col 9 - downstreamES
#          col 10 - downstreamEE
#          col 11 - IJC_SAMPLE_1
#          col 12 - SJC_SAMPLE_2
#

allSorted="all.SE.sorted.txt"
sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $allSE > $allSorted


#
#   step 4 - cut the last two columns from sorted unique union file
#
# cut the last two columns - because they do not represent the sample specific data
allCut="all.SE.sorted.cut.txt"

cut -f 1-10 $allSorted > $allCut

#
#   step 5 - add an ID
#
# add an ID - counting each row
# afterwards we have the following order in the allUnionFile
#
#          col 1 - ID new -- adding now
#          col 2 - GeneID
#          col 3 - geneSymbol
#          col 4 - chr
#          col 5 - strand
#          col 6 - exonStart_0base
#          col 7 - exonEnd
#          col 8 - upstreamES
#          col 9 - upstreamEE
#          col 10 - downstreamES
#          col 11 - downstreamEE

allUnionFile="all.SE.sorted.cut.nl.txt"
nl $allCut > $allUnionFile

#
# step 6 - sort the individual files
#
# sort the individual files (here we keep the sample count)
sortedSEend=".sorted.SE.txt"

for file in $allSEend; do
    name="${file%-*.SE.txt}"
    name_sorted_se=$name$sortedSEend

    echo "file                = " $file
    echo "name                = " $name
    echo "name_sorted_se      = " $name_sorted_se

    sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $file > $name_sorted_se

done

#
# step 7 - normalize the individual files
#
allSortedSE="*.sorted.SE.txt"
normSEend=".sorted.norm.SE.txt"

for file in $allSortedSE; do
    name="${file%-*.sorted.SE.txt}"
    name_norm_se=$name$normSEend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_norm_se        = " $name_norm_se

    awk -f "../bin/match_se.awk" $file $allUnionFile > $name$normSEend

done
