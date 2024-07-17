#!/bin/bash
#---------------------------------------------------------------------
#
# name:  prepareA5SSfiles.sh
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
#   step 1 - extract A5SS coordinages with IJC and SJC
#
#          take the output from rMATS with the coordinates of each
#          splice defined
#          extract the A5SS coordinates IJC and SJC numbers
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
#   step 2 - cat the A5SS files with coordinates and counts together (create a union)
#         
#   step 3 - unique sort master union file (all.A5SS.txt)
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
#   step 9 - create A5SS.coordinates.matrix.txt
#                   A5SS.IJC.w.coordinates.matrix.txt
#                   A5SS.SJC.w.coordinates.matrix.txt
#
#            join all the files together by ID
#            as each file is joined add to the header the name.
#            finish with adding the header to the final matrix
#            all.A5SS.IJC.txt
#            all.A5SS.SJC.txt
#
#   step 10 - create the A5SS.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            
#
# output:  the coordinates and gene identifiers, strand for A5SS
#---------------------------------------------------------------------

cd $1
allA5SS="*.A5SS.MATS.JC.txt"
A5SSend=".A5SS.txt"
tmp="tmp.txt"
PWD=$(pwd)
echo "Current Working Directory is = " $PWD
echo "allA5SS                        = " $allA5SS

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Construct the relative path to the awk script
AWK_SCRIPT="$SCRIPT_DIR/../bin/match_a5ss.awk"

#
#   step 1 - extract coordinates from the A5SS file
#
for file in $allA5SS; do
    name="${file%-*.A5SS.MATS.JC.txt}"
    name_a5ss=$name$A5SSend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_a5ss             = " $name_a5ss
    
    cut -f 2-11,13-14 $file > $name_a5ss

    # remove the header
    tail -n +2 $name_a5ss > $tmp && mv $tmp $name_a5ss
done

#
#   step 2 - create the master union file by cat'ing these files together
#
allA5SSend="*.A5SS.txt"
allA5SS="A5SS.all.txt"

cat $allA5SSend > $allA5SS

#
#   step 3 - unique sort master union file (all.A5SS.txt)
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

allSorted="all.A5SS.sorted.txt"
sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $allA5SS > $allSorted


#
#   step 4 - cut the last two columns from sorted unique union file
#
# cut the last two columns - because they do not represent the sample specific data
allCut="all.A5SS.sorted.cut.txt"

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

a5ssCoordinatesFile="A5SS.coordinates.matrix.txt"
nl $allCut > $a5ssCoordinatesFile

#
# step 6 - sort the individual files
#
# sort the individual files (here we keep the sample count)
sortedA5SSend=".sorted.A5SS.txt"

for file in $allA5SSend; do
    name="${file%.A5SS.txt}"
    name_sorted_a5ss=$name$sortedA5SSend

    echo "file                = " $file
    echo "name                = " $name
    echo "name_sorted_a5ss    = " $name_sorted_a5ss

    sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $file > $name_sorted_a5ss

done

#
# step 7 - normalize the individual files
#
allSortedA5SS="*.sorted.A5SS.txt"
normA5SSend=".sorted.norm.A5SS.txt"

for file in $allSortedA5SS; do
    name="${file%.sorted.A5SS.txt}"
    name_norm_a5ss=$name$normA5SSend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_norm_a5ss      = " $name_norm_a5ss

    # old absolute path
    #    awk -f "/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/match_a5ss.awk" $file $a5ssCoordinatesFile > $name$normA5SSend
    # new relative path
    gawk -f "$AWK_SCRIPT" $file $a5ssCoordinatesFile > $name$normA5SSend
    
done

#   step 8 - break normalized files into IJC, SJC files
#          - coordinate file master will be ID and coordinates as it is
#          - norm IJC file will be ID IJC
#          - norm SJC file will be ID SJC
#
allNormA5SS="*.sorted.norm.A5SS.txt"
ijc=".IJC.txt"
sjc=".SJC.txt"

for file in $allNormA5SS; do
    name="${file%.sorted.norm.A5SS.txt}"
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
#            all.A5SS.IJC.txt
#            all.A5SS.SJC.txt
#
a5ssCoordinatesFile="A5SS.coordinates.matrix.txt"

IJC_matrix="A5SS.IJC.matrix.txt"
IJC_matrix_csv="A5SS.IJC.matrix.csv"
IJC_w_coordinates_matrix="A5SS.IJC.w.coordinates.matrix.txt"
IJC_w_coordinates_matrix_csv="A5SS.IJC.w.coordinates.matrix.csv"
SJC_matrix="A5SS.SJC.matrix.txt"
SJC_matrix_csv="A5SS.SJC.matrix.csv"
SJC_w_coordinates_matrix="A5SS.SJC.w.coordinates.matrix.txt"
SJC_w_coordinates_matrix_csv="A5SS.SJC.w.coordinates.matrix.csv"
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
header_file="A5SS.header.txt"
coordinates_header_file="A5SS.coordinates.header.txt"

header="ID"
coordinates_header="ID	GeneID	geneSymbol	chr	strand	exonStart_0base	exonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE"
tab="	"

#
# Populate the IJC and SJC matrices with the ID from the coordinate file
#
cut -f 1 $a5ssCoordinatesFile > $IJC_matrix
cut -f 1 $a5ssCoordinatesFile > $SJC_matrix

cp $a5ssCoordinatesFile $IJC_w_coordinates_matrix
cp $a5ssCoordinatesFile $SJC_w_coordinates_matrix


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
#   step 10 - create the A5SS.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            

echo "track name=rMATS_A5SS description=\"rMATS A5SS Events DS-AML\"" > A5SS.coordinates.bed

# old absolute path
#     awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_a5ss.awk A5SS.coordinates.matrix.txt >> A5SS.coordinates.bed
# new relative path
gawk -f "$SCRIPT_DIR/../bin/make_bed_a5ss.awk" A5SS.coordinates.matrix.txt >> A5SS.coordinates.bed

