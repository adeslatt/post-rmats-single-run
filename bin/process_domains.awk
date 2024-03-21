#!/usr/bin/awk -f

# Set the field separator to ":"
BEGIN { FS = ":" }

# Initialize the variable for the first time
BEGIN { first_time = 1 }

# Process each line from the protein domain data file
{
    # Assign variables based on the fields from the data file
    protein_name = $1
    domain_name = $2
    aa_position = $3
    domain_sequence = $4

    # Debug: Print the values of the variables
    print "DEBUG: Protein:", protein_name, "Domain Name:", domain_name, "AA Position:", aa_position, "Domain Sequence:", domain_sequence

    # If it's the first line, populate experiment_files array with files in the directory
    if (first_time == 1) {
        # Get the list of experiment files in the specified directory ending with "_linear_aa.fa"
        cmd = "ls -1 " experiment_dir "/*_linear_aa.fa"
        print "DEBUG: Command:", cmd
        while ((cmd | getline experiment_file) > 0) {
	    sub(".*/", "", experiment_file)
	    experiment_files[experiment_file] = 1
            print "DEBUG: Experiment File Added:", experiment_file
        }
        close(cmd)
        first_time = 0
    }

    # Debug: Print the number of experiment files found
    print "DEBUG: Number of Experiment Files:", length(experiment_files)

    # Loop through each experiment file in the directory
    for (experiment_file in experiment_files) {
        # Check if domain_sequence is not empty and experiment_file is not an empty string
        if (domain_sequence != "" && experiment_file != "") {
            # Use grep to search for the domain sequence in the experiment file
            cmd = "grep -c \"" domain_sequence "\" \"" experiment_dir "/" experiment_file "\""
            cmd | getline read_count
            close(cmd)

            # Generate the output filename based on the input experiment file
            output_file = experiment_file "_results.txt"

            # If it's the first time for this experiment file, print the header
            if (!(output_file in output_files)) {
                print "Protein,Domain_Name,AA_position,Domain_Sequence,"experiment_file > output_file
                output_files[output_file] = 1
            }

            # Print the results into the current directory
            print protein_name "," domain_name "," aa_position "," domain_sequence "," read_count >> output_file
            print "DEBUG: Writing Result to Output File:", output_file
            print "DEBUG: Results:", protein_name, domain_name, aa_position, domain_sequence, read_count
        }
    }
}
