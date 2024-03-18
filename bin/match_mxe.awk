# Read file B, each of the new sorted arrays and store the matching columns in an associative array
# an associative array is a key-value pair
# where the key is the unique identifier that yields the value
# in our case the values are the IJC and SJC which are unique to this
# sample
#
# Each individual sorted array is as follows:
#
#          col 1 - GeneID
#          col 2 - geneSymbol
#          col 3 - chr
#          col 4 - strand
#          col 5 - 1stexonStart_0base
#          col 6 - 1stexonEnd
#          col 7 - 2ndexonStart_0base
#          col 8 - 2ndexonEnd
#          col 9 - upstreamES
#          col 10 - upstreamEE
#          col 11 - downstreamES
#          col 12 - downstreamEE
#          col 13 - IJC_SAMPLE_1
#          col 14 - SJC_SAMPLE_1
#

# A is the union file and is as follows
#
#          col 1 - ID
#          col 2 - GeneID
#          col 3 - geneSymbol
#          col 4 - chr
#          col 5 - strand
#          col 6 - 1stexonStart_0base
#          col 7 - 1stexonEnd
#          col 8 - 2ndexonStart_0base
#          col 9 - 2ndexonEnd
#          col 10 - upstreamES
#          col 11 - upstreamEE
#          col 12 - downstreamES
#          col 13 - downstreamEE
#
#
# The match will print out the content of A, the union file, and the Counts from sample B, the individual file of IJC and SJC
#
# this ensures all the individual samples files have the same structure and IDs for making a single matrix
#
NR == FNR {
    OFS = "\t"
  #
  # The key is made up of the columns that uniquely define the SE
  # chr strand 1stexonStart_0base 1stexonEnd 2ndexonStart_0base 2ndexonEnd upstreamES upstreamEE downstreamES downstreamEE
  # col col    col                col        col                col        col        col        col          col
  #  3   4      5                 6          7                  8          9          10         11           12
  #
  # these will become our key
  key = $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12
  #
  # our value is defined as the IJC and SJC which are
  # IJC_SAMPLE_1 SJC_SAMPLE_1
  #  col         col
  #   13          14
  #
  value = $13"\t"$14
  matchArray[key] = value
  next
}

# Check if the key from the union file exists in the individual sample array,
# if yes, output the matching values from B; otherwise, output columns from A with zeros for columns 7 and 8
{
    OFS = "\t"

    key = $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12 OFS $13

    if (key in matchArray) {
	print $0"\t"matchArray[key]
    } else {
	print $0"\t0\t0"
    }
}
