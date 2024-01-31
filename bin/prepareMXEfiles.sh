#!/bin/bash
#---------------------------------------------------------------------
#
# name:  prepareMXEfiles.sh
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
#   step 1 - extract MXE coordinages with IJC and SJC
#
#          take the output from rMATS with the coordinates of each
#          splice defined
#          extract the MXE coordinates IJC and SJC numbers
#          col 1 - ID -- DO NOT NEED/USE
#          col 2 - GeneID
#          col 3 - geneSymbol
#          col 4 - chr
#          col 5 - strand
#          col 6 - 1stExonStart_0base
#          col 7 - 1stExonEnd
#          col 8 - 2ndExonStart_0base
#          col 9 - 2ndExonEnd
#          col 10 - upstreamES
#          col 11 - upstreamEE
#          col 12 - downstreamES
#          col 13 - downstreamEE
#          col 14 - ID - DO NOT NEED/USE
#          col 15 - IJC_SAMPLE_1
#          col 16 - SJC_SAMPLE_2
#
#          remove the header
#
#   step 2 - cat the MXE files with coordinates and counts together (create a union)
#         
#   step 3 - unique sort master union file (all.MXE.txt)
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
#   step 9 - create MXE.coordinates.matrix.txt
#                   MXE.IJC.w.coordinates.matrix.txt
#                   MXE.SJC.w.coordinates.matrix.txt
#
#            join all the files together by ID
#            as each file is joined add to the header the name.
#            finish with adding the header to the final matrix
#            all.MXE.IJC.txt
#            all.MXE.SJC.txt
#
#   step 10 - create the MXE.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            
#
# output:  the coordinates and gene identifiers, strand for MXE
#---------------------------------------------------------------------

cd $1
allMXE="*.MXE.MATS.JC.txt"
MXEend=".MXE.txt"
tmp="tmp.txt"
PWD=$(pwd)
echo "Current Working Directory is = " $PWD
echo "allMXE                        = " $allMXE

#
#   step 1 - extract coordinates from the MXE file
#
for file in $allMXE; do
    name="${file%-*.MXE.MATS.JC.txt}"
    name_mxe=$name$MXEend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_mxe             = " $name_mxe
    
    cut -f 2-11,13-14 $file > $name_mxe

    # remove the header
    tail -n +2 $name_mxe > $tmp && mv $tmp $name_mxe
done

#
#   step 2 - create the master union file by cat'ing these files together
#
allMXEend="*.MXE.txt"
allMXE="MXE.all.txt"

cat $allMXEend > $allMXE

#
#   step 3 - unique sort master union file (all.MXE.txt)
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

allSorted="all.MXE.sorted.txt"
sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $allMXE > $allSorted


#
#   step 4 - cut the last two columns from sorted unique union file
#
# cut the last two columns - because they do not represent the sample specific data
allCut="all.MXE.sorted.cut.txt"

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

mxeCoordinatesFile="MXE.coordinates.matrix.txt"
nl $allCut > $mxeCoordinatesFile

#
# step 6 - sort the individual files
#
# sort the individual files (here we keep the sample count)
sortedMXEend=".sorted.MXE.txt"

for file in $allMXEend; do
    name="${file%.MXE.txt}"
    name_sorted_mxe=$name$sortedMXEend

    echo "file                = " $file
    echo "name                = " $name
    echo "name_sorted_mxe      = " $name_sorted_mxe

    sort -u -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 $file > $name_sorted_mxe

done

#
# step 7 - normalize the individual files
#
allSortedMXE="*.sorted.MXE.txt"
normMXEend=".sorted.norm.MXE.txt"

for file in $allSortedMXE; do
    name="${file%.sorted.MXE.txt}"
    name_norm_mxe=$name$normMXEend
    
    echo "file                = " $file
    echo "name                = " $name
    echo "name_norm_mxe        = " $name_norm_mxe

    awk -f "/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/match_mxe.awk" $file $mxeCoordinatesFile > $name$normMXEend

done

#   step 8 - break normalized files into IJC, SJC files
#          - coordinate file master will be ID and coordinates as it is
#          - norm IJC file will be ID IJC
#          - norm SJC file will be ID SJC
#
allNormMXE="*.sorted.norm.MXE.txt"
ijc=".IJC.txt"
sjc=".SJC.txt"

for file in $allNormMXE; do
    name="${file%.sorted.norm.MXE.txt}"
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
#            all.MXE.IJC.txt
#            all.MXE.SJC.txt
#
mxeCoordinatesFile="MXE.coordinates.matrix.txt"

IJC_matrix="MXE.IJC.matrix.txt"
IJC_matrix_csv="MXE.IJC.matrix.csv"
IJC_w_coordinates_matrix="MXE.IJC.w.coordinates.matrix.txt"
IJC_w_coordinates_matrix_csv="MXE.IJC.w.coordinates.matrix.csv"
SJC_matrix="MXE.SJC.matrix.txt"
SJC_matrix_csv="MXE.SJC.matrix.csv"
SJC_w_coordinates_matrix="MXE.SJC.w.coordinates.matrix.txt"
SJC_w_coordinates_matrix_csv="MXE.SJC.w.coordinates.matrix.csv"
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
header_file="MXE.header.txt"
coordinates_header_file="MXE.coordinates.header.txt"

header="ID"
coordinates_header="ID	GeneID	geneSymbol	chr	strand	1stStart_0base	1stexonEnd	2ndExonStart_0base	2ndExonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE"
tab="	"

#
# Populate the IJC and SJC matrices with the ID from the coordinate file
#
cut -f 1 $mxeCoordinatesFile > $IJC_matrix
cut -f 1 $mxeCoordinatesFile > $SJC_matrix

cp $mxeCoordinatesFile $IJC_w_coordinates_matrix
cp $mxeCoordinatesFile $SJC_w_coordinates_matrix


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
#   step 10 - create the MXE.coordinates.bed file for display
#            in a genome browser, such as the UCSC genome browser
#            

echo "track name=rMATS_MXE description=\"rMATS MXE Events DS-AML\"" > MXE.coordinates.bed
awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_mxe.awk MXE.coordinates.matrix.txt >> MXE.coordinates.bed
