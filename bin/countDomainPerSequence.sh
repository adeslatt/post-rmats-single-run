#!/bin/bash
#
#  program countDomainPerSequence.sh
#
#  purpose:
#   Given open reading frames as linearized (without carriage routine)
#   as amino acid sequences from mRNA experiments (preferabbly long read mRNA)
#   search the protein domains that are included in the measurements.
#
#  why:
#    One set of data is from single cell mRNA long read sequencing.
#    we are checking to see if there are putatively likely alternative splicing
#    isoforms in our mixture.  This would point to potentially different functions
#    and would explain differential expression of the transcript isoforms
#
#
#  
#
#
# Ensure proper usage of the script
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <experiment_dir> <protein_domain_data_file>"
    exit 1
fi

# Assign command line arguments to variables
experiment_dir=$1
data_file=$2

# Call the awk script with experiment directory as an argument
awk -v experiment_dir="$experiment_dir" -f "/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/process_domains.awk" "$data_file"
