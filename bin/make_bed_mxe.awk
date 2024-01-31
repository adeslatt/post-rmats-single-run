#
# make bed mxe awk
#
#  The bed file is constructed in such away that the alternative splicing event
#  may be displayed in a browser such as the UCSC Genome Browser
#
#  The routine follows the elements of style approach, input - process - output precisely as follows:
#
#  INPUT - one of the output files that is generated with the prepareMXEfiles.sh script, MXE.coordinates.matrix.txt
#
#
#     The format of the input is as follows:
#
#     col 1 - ID - unique identifier for the skipped exon event
#     col 2 - GeneID - the ENSG identifier
#     col 3 - geneSymbol - the text word for the gene
#     col 4 - chromsome
#     col 5 - strand
#     col 6 - 1stexonStart_0base - the start base coordinate (zero based) for the exon of interest
#     col 7 - 1stexonEnd - the end base coordinate for the exon of interest
#     col 8 - 2ndexonStart_0base - start of 2nd exon
#     col 9 - 2ndexonEnd
#     col 10 - upstreamES - the start base coordinate (zero based) for the upstream exon
#     col 11 - upstreamEE - the end coordinate for the upstream exon
#     col 12 - downstreamES - the start base coordinate (zero based) for the downstream exon
#     col 13 - downstreamEE - the end coordiante for the downstream exon
#
#  PROCESS -
#    Use the input file to generate the output file.
#    The output file format is described below in the output section
#    processing involves rearrangements:
#
#   INPUT           OUTPUT
#    col 4  becomes col 1 - that is the chromsome from the input becomes the chromsome in the output
#    col 10 becomes col 2 - that is the upstream start exon becomes the thick start
#    col 13 becomes col 3 - the downstream end coordinate becomes the thick end
#    col 1  becomes col 4 - that is the unique identifer for the MXE event is the name
#    0              col 5 - arbitrary score zero
#    col 5  becomes col 6 - the strand
#    col 8  becomes col 7 - repeat of col 2
#    col 11 becomes col 8 - repeat of col 3
#    0 black color  col 9
#    4 exon count   col 10
#    some math (col 11 - col 10),(col 7 - col 6), (col 9 - col 8), (col 11 - col 10), (col 13 - col 12) sizes of the exons upstream, 1st exon, 2ndexon, downstream
#    come more math (col 6 - col 6),(col 8 - col 6),(col 10 - col 6) (col 12 - col 6) - the relative start positions of the exons
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
#    col 2 - is the "thickStart" which is the beginning (this is the upstreamStart coordinate of the upstream exon) of the piece of mRNA we are exploring
#    col 3 - is the "thick end" which is the end of the piece of mRNA we are exploring (the downstreamEnd coordinate of the downstream Exon)
#    col 4 - the name for the piece, in our case we will use the unique identifier generated from the creation of the unified file of all the
#            samples under study
#    col 5 - a score - for us this iz zero
#    col 6 - the strand
#    col 7 - a repeat of col 2
#    col 8 - a repeat of col 3
#    col 9 - color of the display - we set it to black - 0
#    col 10 - the number of exons, for the MXE event this is always 3
#    col 11 - a comma separated list of the lengths of the exons in the following order upstream exon, exon, downstream exon
#    col 12 - a comma separated list of the relative start positions of the exons calculated
#
#  OUTPUT: MXE.coordinates.matrix.bed
#
{
    # trippy for me is that the exons are now rearranged
    # exon2  is one of the two  exons of interest
    # exon3  is the second of two exons of interest
    # exon1 is the upstream exon
    # exon4 is the downstream exon
    exon1 = $11 - $10 #upstream
    exon2 = $7 - $6   #1st exon in MXE 
    exon3 = $9 - $8   #2nd exon in MXE
    exon4 = $13 - $12 #downstream

    # so the transcript start is the upstream start col
    start1 = $10 - $10
    start2 = $6  - $10
    start3 = $8  - $10
    start4 = $12 - $10
    print $4 OFS $10 OFS $13 OFS $1 OFS 0 OFS $5 OFS $10 OFS $13 OFS 0 OFS 4 OFS exon1","exon2","exon3","exon4 OFS start1","start2","start3","start4
}
