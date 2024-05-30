v#
# make bed a3ss awk
#
#  The bed file is constructed in such away that the alternative splicing event
#  may be displayed in a browser such as the UCSC Genome Browser
#
#  The routine follows the elements of style approach, input - process - output precisely as follows:
#
#  INPUT - one of the output files that is generated with the prepareA3SSfiles.sh script, A3SS.coordinates.matrix.txt
#
#
#     The format of the input is as follows:
#
#          col 1 - ID
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
#
#  PROCESS -
#    Use the input file to generate the output file.
#    The output file format is described below in the output section
#    processing involves rearrangements:
#
#   INPUT           OUTPUT
#    col 4   becomes col 1 - that is the chromsome from the input becomes the chromsome in the output
#    col 10  becomes col 2 - that is the flanking exon start becomes the thick start
#    col 9   becomes col 3 - the short exon end coordinate becomes the thick end
#    col 1   becomes col 4 - that is the unique identifer for the A3SS event is the name
#    0               col 5 - arbitrary score zero
#    col 5   becomes col 6 - the strand
#    col 10  becomes col 7 - repeat of col 2
#    col 9   becomes col 8 - repeat of col 3
#    0 black color  col 9
#    3 exon count   col 10
#    some math (col 9 - col 8),(col 7 - col 6), (col 11 - col 10) - size of the exons short exon, long exon, flanking exon
#    come more math (col 6 - col 11),(col 8 - col 11),(col 10 - col 11) - the relative start positions of the exons
#
#
#  OUTPUT - 
#    The format of the file that we will produces follows the Bed File format
#    Track name may be obtained from command line input - likely programmatically in the future
#    The description also similiarly - default will be Gene Name and splicing type.
#
#    Example is below: (generated with help from chatGPT-3.5
#      track name=my_coding_mRNA description="My coding mRNA with multiple exons"
#      chr6    36596760    36598983    MYGENE    0    +    36596760    36598983    0    3    208,98,135    0,1668,2088
#
#    col 1 - the chromsome id
#    col 2 - is the "thickStart" which is the beginning (this is the short exon (Start coordinate of the short exon) of the piece A3SS
#    col 3 - is the "thick end" which is the end of the piece of mRNA we are exploring (the flanking exon start coordinate of the flanking Exon)
#    col 4 - the name for the piece, in our case we will use the unique identifier generated from the creation of the unified file of all the
#            samples under study
#    col 5 - a score - for us this iz zero
#    col 6 - the strand
#    col 7 - a repeat of col 2
#    col 8 - a repeat of col 3
#    col 9 - color of the display - we set it to black - 0
#    col 10 - the number of exons, for the A3SS event this is always 3
#    col 11 - a comma separated list of the lengths of the exons in the following order short exon, long exon, flanking exon
#    col 12 - a comma separated list of the relative start positions of the exons calculated
#
#  OUTPUT: A3SS.coordinates.bed
#

BEGIN {

    # make sure the default delimiter for output printing is a tab
    OFS = "\t"
}

# Assuming the header line is being read but not printed here
NR == 1 {
    next;  # Skip processing and move to the next line
}

{
    # with A3SS - the short exon, refers to the A3ss making the 3 prime exon in a junction smaller, thereby making it a greater distance in fact
    #    then the longer 
    shortExon    = $9 - $8
    longExon     = $7 - $6
    flankingExon = $11 - $10

    # strand
    strand       = $5

    # so the gene symbol that gets printed out does not have quotes on it - lets strip them
    gsub(/"/, "", $3)

    # so the transcript start is the flanking exon start col
    # we also have to order the exons on the chromosome appropriately
    # on the positive strand
    #  flanking is first, long is second and short is third
    if (strand == "+") {
	flankingExonStart = $10 - $10
	longExonStart     = $6 - $10
	shortExonStart    = $8 - $10

	print $4 OFS $10 OFS $9 OFS $3"_"$1 OFS 0 OFS $5 OFS $10 OFS $9 OFS 0 OFS 3 OFS flankingExon","longExon","shortExon OFS flankingExonStart","longExonStart","shortExonStart
    } else {
    # on the negative strand
    #  flanking exon is the highest distance - e.g. the exact opposite since we are going in reverse
    #  short Exon is first, long is second and flanking is third
	shortExonStart    = $8 - $8
	longExonStart     = $6 - $8
	flankingExonStart = $10 - $8

	print $4 OFS $8 OFS $11 OFS $3"_"$1 OFS 0 OFS $5 OFS $8 OFS $11 OFS 0 OFS 3 OFS shortExon","longExon","flankingExon OFS shortExonStart","longExonStart","flankingExonStart
    }
}
