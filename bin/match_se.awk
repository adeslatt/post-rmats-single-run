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
#          col 5 - exonStart_0base
#          col 6 - exonEnd
#          col 7 - upstreamES
#          col 8 - upstreamEE
#          col 9 - downstreamES
#          col 10 - downstreamEE
#          col 11 - IJC_SAMPLE_1
#          col 12 - SJC_SAMPLE_1
#

# A is the union file and is as follows
#
#          col 1 - ID
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
#
#
# The match will print out the content of A, the union file, and the Counts from sample B, the individual file of IJC and SJC
#
# this ensures all the individual samples files have the same structure and IDs for making a single matrix
#
NR == FNR {
  #
  # The key is made up of the columns that uniquely define the SE
  # chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE
  # col col    col             col     col        col        col          col
  #  3   4      5               6       7          8          9           10
  #
  # these will become our key
  key = $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10
  #
  # our value is defined as the IJC and SJC which are
  # IJC_SAMPLE_1 SJC_SAMPLE_1
  #  col         col
  #   11          12
  #
  value = $11"\t"$12
  matchArray[key] = value
  next
}

# Check if the key from the union file exists in the individual sample array,
# if yes, output the matching values from B; otherwise, output columns from A with zeros for columns 7 and 8
{
  key = $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11
  if (key in matchArray) {
    print $0"\t"matchArray[key]
  } else {
    print $0"\t0\t0"
  }
}
