#!/bin/bash
#---------------------------------------------------------------------
#
# name:  prepareRIfiles.sh
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
#   step 1 - extract RI coordinages with IJC and SJC
#
#          take the output from rMATS with the coordinates of each
#          splice defined
#          extract the RI coordinates IJC and SJC numbers
#          col 1 - ID -- DO NOT NEED/USE
#          col 2 - GeneID
#          col 3 - geneSymbol
#          col 4 - chr
#          col 5 - strand
#          col 6 - riexonStart_0base
#          col 7 - riexonEnd
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
#   step 2 - cat the RI files with coordinates and counts together (create a union)
#         
#   step 3 - unique sort master union file (all.RI.txt)
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
#   step 9 - create RI.coordinates.matrix.txt
#                   RI.IJC.w.coordinates.matrix.txt
#                   RI.SJC.w.coordinates.matrix.txt
#
#            join all the files together by ID
#            as each file is joined add to the header the name.
#            finish with adding the header to the final matrix
#            all.RI.IJC.txt
#            all.RI.SJC.txt
#
#   step 10 - create the RI.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            
#
# output:  the coordinates and gene identifiers, strand for RI
#---------------------------------------------------------------------

cd $1
allRI="*.RI.MATS.JC.txt"
RIend=".RI.txt"
tmp="tmp.txt"
PWD=$(pwd)
echo "Current Working Directory is = " $PWD
echo "allRI                        = " $allRI

#
#   step 1 - extract coordinates from the RI file
#
for file in $allRI; do
    name="${file%-*.RI.MATS.JC.txt}"
    name_ri=$name$RIend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_ri             = " $name_ri

    # col 1,12 are the ids which we won't use
    # col 13,14 are the IJC and SJC which we will use
    cut -f 2-11,13-14 $file > $name_ri

    # remove the header
    tail -n +2 $name_ri > $tmp && mv $tmp $name_ri
done

#
#   step 2 - create the master union file by cat'ing these files together
#
allRIend="*.RI.txt"
allRI="RI.all.txt"

cat $allRIend > $allRI

#
#   step 3 - unique sort master union file (all.RI.txt)
#
# sort unique
#   our column structure has changed
#          col 1 - GeneID
#          col 2 - geneSymbol
#          col 3 - chr
#          col 4 - strand
#          col 5 - riexonStart_0base
#          col 6 - riexonEnd
#          col 7 - upstreamES
#          col 8 - upstreamEE
#          col 9 - downstreamES
#          col 10 - downstreamEE
#          col 11 - IJC_SAMPLE_1
#          col 12 - SJC_SAMPLE_2
#

allSorted="all.RI.sorted.txt"

# sort starting with col 3 - chromosome
#                    col 4 - strand
#                    col 5 - riexonStart_0base
#                    col 6 - riexonEnd
#                    col 7 - upstreamES
#                    col 8 - upstreamEE
#                    col 9 - downstreamES
#                    col 10 - downstreamEE
#
sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $allRI > $allSorted


#
#   step 4 - cut the last two columns from sorted unique union file
#
# cut the last two columns - because they do not represent the sample specific data
allCut="all.RI.sorted.cut.txt"

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
#          col 6 - riexonStart_0base
#          col 7 - riexonEnd
#          col 8 - upstreamES
#          col 9 - upstreamEE
#          col 10 - downstreamES
#          col 11 - downstreamEE

riCoordinatesFile="RI.coordinates.matrix.txt"

# the function nl will add an id to each of the sorted lines - our new master union file
nl $allCut > $riCoordinatesFile

#
# step 6 - sort the individual files
#
# sort the individual files (here we keep the sample count)
sortedRIend=".sorted.RI.txt"

for file in $allRIend; do
    name="${file%.RI.txt}"
    name_sorted_ri=$name$sortedRIend

    echo "file                = " $file
    echo "name                = " $name
    echo "name_sorted_ri      = " $name_sorted_ri

    sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $file > $name_sorted_ri

done

#
# step 7 - normalize the individual files
#
allSortedRI="*.sorted.RI.txt"
normRIend=".sorted.norm.RI.txt"

for file in $allSortedRI; do
    name="${file%.sorted.RI.txt}"
    name_norm_ri=$name$normRIend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_norm_ri        = " $name_norm_ri

    awk -f "/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/match_ri.awk" $file $riCoordinatesFile > $name$normRIend

done

#   step 8 - break normalized files into IJC, SJC files
#          - coordinate file master will be ID and coordinates as it is
#          - norm IJC file will be ID IJC
#          - norm SJC file will be ID SJC
#
allNormRI="*.sorted.norm.RI.txt"
ijc=".RI.IJC.txt"
sjc=".RI.SJC.txt"

for file in $allNormRI; do
    name="${file%.sorted.norm.RI.txt}"
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
#            all.RI.IJC.txt
#            all.RI.SJC.txt
#
riCoordinatesFile="RI.coordinates.matrix.txt"

IJC_matrix="RI.IJC.matrix.txt"
IJC_matrix_csv="RI.IJC.matrix.csv"
IJC_w_coordinates_matrix="RI.IJC.w.coordinates.matrix.txt"
IJC_w_coordinates_matrix_csv="RI.IJC.w.coordinates.matrix.csv"
SJC_matrix="RI.SJC.matrix.txt"
SJC_matrix_csv="RI.SJC.matrix.csv"
SJC_w_coordinates_matrix="RI.SJC.w.coordinates.matrix.txt"
SJC_w_coordinates_matrix_csv="RI.SJC.w.coordinates.matrix.csv"
allIJC="*.RI.IJC.txt"
allSJC="*.RI.SJC.txt"
IJCend=".RI.IJC.txt"
SJCend=".RI.SJC.txt"

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
header_file="RI.header.txt"
coordinates_header_file="RI.coordinates.header.txt"

header="ID"
coordinates_header="ID	GeneID	geneSymbol	chr	strand	riexonStart_0base	riexonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE"
tab="	"

#
# Populate the IJC and SJC matrices with the ID from the coordinate file
#
cut -f 1 $riCoordinatesFile > $IJC_matrix
cut -f 1 $riCoordinatesFile > $SJC_matrix

cp $riCoordinatesFile $IJC_w_coordinates_matrix
cp $riCoordinatesFile $SJC_w_coordinates_matrix


for file in $allIJC; do
    name="${file%.RI.IJC.txt}"
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
#   step 10 - create the RI.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            

echo "track name=rMATS_RI description=\"rMATS RI Events DS-AML\"" > RI.coordinates.bed
awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_ri.awk RI.coordinates.matrix.txt >> RI.coordinates.bed
