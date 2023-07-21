# Read file B and store the matching columns in an associative array
# an associative array is a key-value pair
# where the key is the unique identifier that yields the value
# in our case the values are the IJC and SJC which are unique to this
# sample
NR == FNR {
  # the key is made up of the columns that uniquely define the SE
  # chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE
  # col col    col             col     col        col        col          col
  #  4   5      6               7       8          9          10           11
  # these will become our key
  key = $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11
  # our value is defined as the IJC and SJC which are
  # IJC_SAMPLE_1 SJC_SAMPLE_1
  #  col         col
  #   13          14
  #
  value = $13 OFS $14
  matchArray[key] = value
  next
}

# Check if the key exists in the array, if yes, output the matching values from B; otherwise, output columns from A with zeros for columns 7 and 8
{
  key = $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11
  if (key in matchArray) {
    print $0, matchArray[key]
  } else {
    print $0, "0 0"
  }
}
