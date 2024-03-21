#!/bin/bash

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
