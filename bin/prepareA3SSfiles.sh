#!/bin/bash
#---------------------------------------------------------------------
#
# name:  prepareA3SSfiles.sh
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
#   step 1 - extract A3SS coordinages with IJC and SJC
#
#          take the output from rMATS with the coordinates of each
#          splice defined
#          extract the A3SS coordinates IJC and SJC numbers
#
#          col 1 - ID -- DO NOT NEED/USE
#          col 2 - GeneID
#          col 3 - geneSymbol
#          col 4 - chr
#          col 5 - strand
#          col 6 - longExonStart_0base
#          col 7 - longExonEnd
#          col 8 - shortES
#          col 9 - shortEE
#          col 10 - flankingES
#          col 11 - flankingEE
#          col 12 - ID - DO NOT NEED/USE
#          col 13 - IJC_SAMPLE_1
#          col 14 - SJC_SAMPLE_2
#
#          remove the header
#
#   step 2 - cat the A3SS files with coordinates and counts together (create a union)
#         
#   step 3 - unique sort master union file (all.A3SS.txt)
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
#   step 9 - create A3SS.coordinates.matrix.txt
#                   A3SS.IJC.w.coordinates.matrix.txt
#                   A3SS.SJC.w.coordinates.matrix.txt
#
#            join all the files together by ID
#            as each file is joined add to the header the name.
#            finish with adding the header to the final matrix
#            all.A3SS.IJC.txt
#            all.A3SS.SJC.txt
#
#   step 10 - create the A3SS.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            
#
# output:  the coordinates and gene identifiers, strand for A3SS
#---------------------------------------------------------------------

cd $1
allA3SS="*.A3SS.MATS.JC.txt"
A3SSend=".A3SS.txt"
tmp="tmp.txt"
PWD=$(pwd)
echo "Current Working Directory is = " $PWD
echo "allA3SS                        = " $allA3SS

#
#   step 1 - extract coordinates from the A3SS file
#
for file in $allA3SS; do
    name="${file%-*.A3SS.MATS.JC.txt}"
    name_a3ss=$name$A3SSend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_a3ss             = " $name_a3ss
    
    cut -f 2-11,13-14 $file > $name_a3ss

    # remove the header
    tail -n +2 $name_a3ss > $tmp && mv $tmp $name_a3ss
done

#
#   step 2 - create the master union file by cat'ing these files together
#
allA3SSend="*.A3SS.txt"
allA3SS="A3SS.all.txt"

cat $allA3SSend > $allA3SS

#
#   step 3 - unique sort master union file (all.A3SS.txt)
#
# sort unique
#   our column structure has changed
#          col 1 - GeneID
#          col 2 - geneSymbol
#          col 3 - chr
#          col 4 - strand
#          col 5 - longExonStart_0base
#          col 6 - longExonEnd
#          col 7 - shortES
#          col 8 - shortEE
#          col 9 - flankingES
#          col 10 - flankingEE
#          col 11 - IJC_SAMPLE_1
#          col 12 - SJC_SAMPLE_2

allSorted="all.A3SS.sorted.txt"
sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $allA3SS > $allSorted


#
#   step 4 - cut the last two columns from sorted unique union file
#
# cut the last two columns - because they do not represent the sample specific data
allCut="all.A3SS.sorted.cut.txt"

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

a3ssCoordinatesFile="A3SS.coordinates.matrix.txt"
nl $allCut > $a3ssCoordinatesFile

#
# step 6 - sort the individual files
#
# sort the individual files (here we keep the sample count)
sortedA3SSend=".sorted.A3SS.txt"

for file in $allA3SSend; do
    name="${file%.A3SS.txt}"
    name_sorted_a3ss=$name$sortedA3SSend

    echo "file                = " $file
    echo "name                = " $name
    echo "name_sorted_a3ss    = " $name_sorted_a3ss

    sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $file > $name_sorted_a3ss

done

#
# step 7 - normalize the individual files
#
allSortedA3SS="*.sorted.A3SS.txt"
normA3SSend=".sorted.norm.A3SS.txt"

for file in $allSortedA3SS; do
    name="${file%.sorted.A3SS.txt}"
    name_norm_a3ss=$name$normA3SSend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_norm_a3ss      = " $name_norm_a3ss

    awk -f "/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/match_a3ss.awk" $file $a3ssCoordinatesFile > $name$normA3SSend

done

#   step 8 - break normalized files into IJC, SJC files
#          - coordinate file master will be ID and coordinates as it is
#          - norm IJC file will be ID IJC
#          - norm SJC file will be ID SJC
#
allNormA3SS="*.sorted.norm.A3SS.txt"
ijc=".IJC.txt"
sjc=".SJC.txt"

for file in $allNormA3SS; do
    name="${file%.sorted.norm.A3SS.txt}"
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
#            all.A3SS.IJC.txt
#            all.A3SS.SJC.txt
#
a3ssCoordinatesFile="A3SS.coordinates.matrix.txt"

IJC_matrix="A3SS.IJC.matrix.txt"
IJC_matrix_csv="A3SS.IJC.matrix.csv"
IJC_w_coordinates_matrix="A3SS.IJC.w.coordinates.matrix.txt"
IJC_w_coordinates_matrix_csv="A3SS.IJC.w.coordinates.matrix.csv"
SJC_matrix="A3SS.SJC.matrix.txt"
SJC_matrix_csv="A3SS.SJC.matrix.csv"
SJC_w_coordinates_matrix="A3SS.SJC.w.coordinates.matrix.txt"
SJC_w_coordinates_matrix_csv="A3SS.SJC.w.coordinates.matrix.csv"
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
header_file="A3SS.header.txt"
coordinates_header_file="A3SS.coordinates.header.txt"

header="ID"
coordinates_header="ID	GeneID	geneSymbol	chr	strand	exonStart_0base	exonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE"
tab="	"

#
# Populate the IJC and SJC matrices with the ID from the coordinate file
#
cut -f 1 $a3ssCoordinatesFile > $IJC_matrix
cut -f 1 $a3ssCoordinatesFile > $SJC_matrix

cp $a3ssCoordinatesFile $IJC_w_coordinates_matrix
cp $a3ssCoordinatesFile $SJC_w_coordinates_matrix


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
# add a csv version for files
#
sed  's/ /,/g' < $IJC_matrix > $IJC_matrix_csv
sed  's/ /,/g' < $SJC_matrix > $SJC_matrix_csv

sed 's/ /,/g' < $IJC_w_coordinates_matrix > $IJC_w_coordinates_matrix_csv
sed 's/ /,/g' < $SJC_w_coordinates_matrix > $SJC_w_coordinates_matrix_csv
#
#   step 10 - create the A3SS.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            

echo "track name=rMATS_A3SS description=\"rMATS A3SS Events DS-AML\"" > A3SS.coordinates.bed
awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_a3ss.awk A3SS.coordinates.matrix.txt >> A3SS.coordinates.bed
