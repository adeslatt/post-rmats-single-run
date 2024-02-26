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
# Define the relative directory where your experimental files are located
#experiment_dir="/path/to/experiment_directory"
experiment_dir=$1
# Define the relative location file containing protein and domain data
#  (one pair per line, in the format "ProteinName:DomainSequence")
#  FUTURE THOUGHTS: Replace this with an API.
#

#data_file="/path/to/protein_domain_data.txt"
data_file=$2

# Loop through each line in the data file
first_time=1
while IFS= read -r line; do
    # Split the line into protein name and domain sequence using ":" as a delimiter
    IFS=":" read -r protein_name domain_name aa_position domain_sequence <<< "$line"
    
    # Loop through each experiment file in the directory
    for experiment_file in "$experiment_dir"/*_linear_aa.fa; do
        # Use grep to search for the domain sequence in the experiment file
        read_count=$(grep -c "$domain_sequence" "$experiment_file")
        
        # Generate the output filename based on the input experiment file
        output_file="${experiment_file##*/}_results.txt"
        if (($first_time == 1));  then
	    echo "Protein Domain_Name AA_position Domain_Sequence $experiment_file" > "$output_file"
	    first_time=0
	fi
        echo "$protein_name $domain_name $aa_position $domain_sequence $read_count" >> $output_file
    done
done < "$data_file"
