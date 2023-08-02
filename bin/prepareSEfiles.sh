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
#   step 8 - break normalized files into IJC, SJC files
#          - coordinate file master will be ID and coordinates as it is
#          - norm IJC file will be ID IJC
#          - norm SJC file will be ID SJC
#
#   step 9 - create SE.coordinates.matrix.txt
#                   SE.IJC.w.coordinates.matrix.txt
#                   SE.SJC.w.coordinates.matrix.txt
#
#            join all the files together by ID
#            as each file is joined add to the header the name.
#            finish with adding the header to the final matrix
#            all.SE.IJC.txt
#            all.SE.SJC.txt
#
#   step 10 - create the SE.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            
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

seCoordinatesFile="SE.coordinates.matrix.txt"
nl $allCut > $seCoordinatesFile

#
# step 6 - sort the individual files
#
# sort the individual files (here we keep the sample count)
sortedSEend=".sorted.SE.txt"

for file in $allSEend; do
    name="${file%.SE.txt}"
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
    name="${file%.sorted.SE.txt}"
    name_norm_se=$name$normSEend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_norm_se        = " $name_norm_se

    awk -f "../bin/match_se.awk" $file $seCoordinatesFile > $name$normSEend

done

#   step 8 - break normalized files into IJC, SJC files
#          - coordinate file master will be ID and coordinates as it is
#          - norm IJC file will be ID IJC
#          - norm SJC file will be ID SJC
#
allNormSE="*.sorted.norm.SE.txt"
ijc=".IJC.txt"
sjc=".SJC.txt"

for file in $allNormSE; do
    name="${file%.sorted.norm.SE.txt}"
    name_ijc=$name$ijc
    name_sjc=$name$sjc

    echo "file                = " $file
    echo "name                = " $name
    echo "name_ijc            = " $name_ijc
    echo "name_sjc            = " $name_sjc

    cut -f 1,12 $file > $name_ijc
    cut -f 1,13 $file > $name_sjc

done

    
#   step 9 - create final matrix
#            join all the files together by ID
#            as each file is joined add to the header the name.
#            finish with adding the header to the final matrix
#            all.SE.IJC.txt
#            all.SE.SJC.txt
#
seCoordinatesFile="SE.coordinates.matrix.txt"

IJC_matrix="SE.IJC.matrix.txt"
IJC_w_coordinates_matrix="SE.IJC.w.coordinates.matrix.txt"
SJC_matrix="SE.SJC.matrix.txt"
SJC_w_coordinates_matrix="SE.SJC.w.coordinates.matrix.txt"
allIJC="*.IJC.txt"
allSJC="*.SJC.txt"
IJCend=".IJC.txt"
SJCend=".SJC.txt"

#
# need tmp files for temporality
#
tmp_IJC="tmp_IJC.txt"
tmp_SJC="tmp_SJC.txt"
tmp_coordinates_IJC="tmp_coord_IJC.txt"
tmp_coordinates_SJC="tmp_coord_SJC.txt"


#
# headers
#
header_file="SE.header.txt"
coordinates_header_file="SE.coordinates.header.txt"

header="ID"
coordinates_header="ID	GeneID	geneSymbol	chr	strand	exonStart_0base	exonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE"
tab="	"

#
# Populate the IJC and SJC matrices with the ID from the coordinate file
#
cut -f 1 $seCoordinatesFile > $IJC_matrix
cut -f 1 $seCoordinatesFile > $SJC_matrix

cp $seCoordinatesFile $IJC_w_coordinates_matrix
cp $seCoordinatesFile $SJC_w_coordinates_matrix


for file in $allIJC; do
    name="${file%.IJC.txt}"
    header=$header$tab$name
    coordinates_header=$coordinates_header$tab$name

    IJCfile=$name$IJCend
    SJCfile=$name$SJCend
    
    echo "IJCfile             = " $IJCfile
    echo "SJCfile             = " $SJCfile
    echo "name                = " $name
    echo "header              = " $header
    echo "coordinates_header  = " $coordinates_header
    
    join -1 1 -2 1 $IJC_matrix $IJCfile > $tmp_IJC
    join -1 1 -2 1 $SJC_matrix $SJCfile > $tmp_SJC

    join -1 1 -2 1 $IJC_w_coordinates_matrix $IJCfile > $tmp_coordinates_IJC
    join -1 1 -2 1 $SJC_w_coordinates_matrix $SJCfile > $tmp_coordinates_SJC
    
    cp $tmp_IJC $IJC_matrix
    cp $tmp_SJC $SJC_matrix

    cp $tmp_coordinates_IJC $IJC_w_coordinates_matrix
    cp $tmp_coordinates_SJC $SJC_w_coordinates_matrix
    
done

#
# Add header
#
echo $header > $header_file
echo $coordinates_header > $coordinates_header_file

cat $header_file $IJC_matrix > $tmp_IJC
cat $header_file $SJC_matrix > $tmp_SJC

cat $coordinates_header_file $IJC_w_coordinates_matrix > $tmp_coordinates_IJC
cat $coordinates_header_file $SJC_w_coordinates_matrix > $tmp_coordinates_SJC

cp $tmp_IJC $IJC_matrix
cp $tmp_SJC $SJC_matrix

cp $tmp_coordinates_IJC $IJC_w_coordinates_matrix
cp $tmp_coordinates_SJC $SJC_w_coordinates_matrix

#
# clean up
#
rm $tmp_IJC
rm $tmp_SJC

rm $tmp_coordinates_IJC
rm $tmp_coordinates_SJC

#
#   step 10 - create the SE.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            

echo "track name=rMATS_SE description=\"rMATS SE Events DS-AML\"" > SE.coordinates.bed
awk -f ../bin/make_bed_se.awk SE.coordinates.matrix.txt >> SE.coordinates.bed
