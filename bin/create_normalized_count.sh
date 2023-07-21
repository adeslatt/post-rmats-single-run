#!/bin/bash
#---------------------------------------------------------------------
#
# name:  create_normalized_SE_count.sh
# purpose: take union of all SE coordinates and now match to coordinates in sample specific count.
# input:   $1 = current directory
#          $2 = name of the file containing the union of SE coordinates with unique ID
#
# Process: loop through the union of all SE file
#               union_se_line_coordinates there are 6 of them -- the exon in question, and the exon upstream and the exon downstream
#          loop through the files that are the contributors to the union
#               se_line_coordinates - is there a match.
#               if there is a match write out all_se_line_coordinate with IJC and SJC counts
#               else write out union_se_coordinates with 0 for IJC and 0 for SJC
#
# output:  matched coordinates and ID with IJC and SJC for Sample
#          Output two matrices, IJC and SJC, with Column Name equal to the Subject_Sample
#  
#---------------------------------------------------------------------

# change to the current directory
PWD=$1
cd $PWD

# get the union SE file
unionSE=$2

# what file to compare
specificSE=$3

#
# Global Variables
#
#allSE="*.sorted.SE.txt"

normSEend=".norm.SE.txt"
tab="	"
nullIJC=0
nullSJC=0

header="ID"$tab"GeneID"$tab"geneSymbol"$tab"chr"$tab"strand"$tab"exonStart_0base"$tab"exonEnd"$tab"upstreamES"$tab"upstreamEE"$tab"downstreamES"$tab"downstreamEE"$tab"IJC"$tab"SJC"


echo "Current Working Directory is = " $PWD
echo "unionSE file is              = " $unionSE
echo "specificSE                   = " $specificSE

#
# two loops - outer loop is the unionSE file
# inner loop is each of the individual SE files.
#
# note that to use the cut command - you need to use 'verbatim' or a quoted insert
# which is accomplished with Ctrl + V and then Tab
# loop through -- in emacs this is done as Ctrl + Q and then tab
#
first=1
count=0
SEcount=0
outSEcount=0
zero=0
one=1
name="${specificSE%.sorted.SE.txt}"
normSE=$name$normSEend

echo "SE name           = " $name
echo "SE norm file name = " $normSE
#echo "count is " $count

#
# First read a line from the union of all the SE files
#
while read -r line
do
    #
    # Skip the header
    #
    if [ "$count" -eq "$zero" ];
    then
	read -r line
	count=1
    fi
    
    # set the marker that match is found
    match=0

    ID="$(cut -d"	" -f1 <<<"$line")"
    GeneID="$(cut -d"	" -f2 <<<"$line")"
    geneSymbol="$(cut -d"	" -f3 <<<"$line")"
    chr="$(cut -d"	" -f4 <<<"$line")"
    strand="$(cut -d"	" -f5 <<<"$line")"
    exonStart_0base="$(cut -d"	" -f6 <<<"$line")"
    exonEnd="$(cut -d"	" -f7 <<<"$line")"
    upstreamES="$(cut -d"	" -f8 <<<"$line")"
    upstreamEE="$(cut -d"	" -f9 <<<"$line")"
    downstreamES="$(cut -d"	" -f10 <<<"$line")"
    downstreamEE="$(cut -d"	" -f11 <<<"$line")"

    
#    echo "Line to match is " $line
#
# Debug code to ensure proper parsing
#
#    echo "ID              = " $ID
#    echo "GeneID          = " $GeneID
#    echo "geneSymbol      = " $geneSymbol
#    echo "chr             = " $chr
#    echo "strand          = " $strand
#    echo "exonStart_0base = " $exonStart_0base
#    echo "exonEnd         = " $exonEnd
#    echo "upstreamES      = " $upstreamES
#    echo "upstreamEE      = " $upstreamES
#    echo "downstreamES    = " $downstreamES
#    echo "downstreamEE    = " $downstreamEE

    #
    # now compare with each line of each specific SE file
    # note that the ID of the specific SE file does not relate to our unionSE and we will not use it
    # so it is not read in and is skipped
    #
    # first we need to create the normalized output file to hold the expanded coordinates
    
    SEcount=0
    while read -r SEline
    do
	
    	#
    	# Skipped the header
    	#
        if [ "$SEcount" -eq "$zero" ];
        then
    	   	 read -r SEline
    	   	 SEcount=1
        fi
   
    	SEGeneID="$(cut -d"	" -f2 <<<"$SEline")"
    	SEgeneSymbol="$(cut -d"	" -f3 <<<"$SEline")"
    	SEchr="$(cut -d"	" -f4 <<<"$SEline")"
    	SEstrand="$(cut -d"	" -f5 <<<"$SEline")"
    	SEexonStart_0base="$(cut -d"	" -f6 <<<"$SEline")"
    	SEexonEnd="$(cut -d"	" -f7 <<<"$SEline")"
    	SEupstreamES="$(cut -d"	" -f8 <<<"$SEline")"
    	SEupstreamEE="$(cut -d"	" -f9 <<<"$SEline")"
    	SEdownstreamES="$(cut -d"	" -f10 <<<"$SEline")"
    	SEdownstreamEE="$(cut -d"	" -f11 <<<"$SEline")"
    	SEIJC="$(cut -d"	" -f13 <<<"$SEline")"
    	SESJC="$(cut -d"	" -f14 <<<"$SEline")"

#    	 echo "SEline is " $SEline   	
#   	 echo "SEGeneID          = " $SEGeneID
#   	 echo "SEgeneSymbol      = " $SEgeneSymbol
#   	 echo "SEchr             = " $SEchr
#   	 echo "SEstrand          = " $SEstrand
#   	 echo "SEexonStart_0base = " $SEexonStart_0base
#   	 echo "SEexonEnd         = " $SEexonEnd
#   	 echo "SEupstreamES      = " $SEupstreamES
#   	 echo "SEupstreamEE      = " $SEupstreamES
#   	 echo "SEdownstreamES    = " $SEdownstreamES
#   	 echo "SEdownstreamEE    = " $SEdownstreamEE
#   	 echo "SEIJC             = " $SEIJC
#   	 echo "SESJC             = " $SESJC

    	#
        # be sure to add spaces after '[' or bash interprets as a command the variable
        #
    	# Do we have a match to the union SE file and our specific SE file ?
    	#
    	# if so break out of the loop and write out the line.
    	#
	if [[ "$exonStart_0base" -eq "$SEexonStart_0base" ]];
	then
	    if [[ "$exonEnd"         -eq "$SEexonEnd"         ]] &&
    	   	    [[ "$upstreamES"      -eq "$SEupstreamES"      ]] &&
    	   	    [[ "$upstreamEE"      -eq "$SEupstreamEE"      ]] &&
    	   	    [[ "$downstreamES"    -eq "$SEdownstreamES"    ]] &&
    	   	    [[ "$downstreamEE"    -eq "$SEdownstreamEE"    ]];
	    then
		match=1
#		echo "we have a match $SEline = " $SEline
#		echo "                   line = " $line
    		break
	    fi
	else
	    if [[ "$exonStart_0base" -gt "$SEexonStart_0base" ]];
	    then
		match=0
#		echo "no match lets print"
		break
	    fi
	fi

    done < $specificSE

    
    # be sure to add spaces after '[' or bash interprets as a command the variable
    #
    if [[ "$match" -eq "$one" ]];
    then
    	normSEline=$ID$tab$GeneID$tab$geneSymbol$tab$chr$tab$strand$tab$exonStart_0base$tab$exonEnd$tab$upstreamES$tab$upstreamEE$tab$downstreamES$tab$downstreamEE$tab$SEIJC$tab$SESJC
	match=0
    else
    	normSEline=$ID$tab$GeneID$tab$geneSymbol$tab$chr$tab$strand$tab$exonStart_0base$tab$exonEnd$tab$upstreamES$tab$upstreamEE$tab$downstreamES$tab$downstreamEE$tab$nullIJC$tab$nullSJC
    fi

    if [ "$outSEcount" -eq "$zero" ];
    then
    	echo $header  > $normSE
	echo $normSEline >> $normSE
	outSEcount=1
    else
	echo $normSEline >> $normSE
    fi
	     
done < $unionSE

