#!/usr/bin/env python

import argparse
import sys
import os

# results is global and actions on the dictionary universal

def process_file(file_path, results):
    print(f"process_file: file_path is {file_path}")

    with open(file_path, "r") as file:
        header = file.readline().strip().split(' ')
        # get the last element in the header
        experiment_name = file_path
        

        for line in file:
            data = line.strip().split(' ')
            
            # Here, data[:4] is a list slice that extracts the first four elements (from index 0 to 3)
            key = tuple(data[:4]) + (experiment_name,)

            # the main ocmponent in the matrix comes from the fifth element aka data[4]
            value = int(data[4])
            print(f"process_file: data is {data[4]}")

            if key not in results:
                results[key] = {}

            # Update results with the values
            results[key] = value
            print(f"process_file: results started with key {key} is {results[key]}")

def main(input_directory, output_matrix):
    print(f"main:Input Directory: {input_directory}")
    print(f"main:Output Matrix: {output_matrix}")


    files = [f for f in os.listdir(input_directory) if f.endswith(".txt")]
    # Sort files for consistent order                                                                                                                      
    files.sort()

    results = {}

    print(f"main: number of files is {len(files)}")
    
    for file in files:
        
        file_path = os.path.join(input_directory, file)
        process_file(file_path, results)
        print(f"main: file is {file_path}")
        print(f"main: results data dictionary (dd) is {results}")

    with open(output_matrix, 'w') as output_file:
        # Writing header                                                                                                                                   
        header = ",".join(["Protein", "Domain_Name", "AA_position", "Domain_Sequence"] + files)
        output_file.write(header + "\n")
        # Writing data

        key_list = []
        
        for key, value in results.items():
            print(f"main: key   is {key}")
            print(f"main: value is {value}")

            # now initiate the row with the tuple
            key_defining_columns = key[:4]

            if key_defining_columns not in key_list:
                key_list.append(key_defining_columns)
                row = ",".join(map(str, key_defining_columns))
                print(f"main: row tuple is {row}")
                for file in files:
                    file_path = os.path.join(input_directory, file)
                    made_key = tuple(key[:4]) + (file_path,)
                    print(f"main: made_key is {made_key}")
                    made_value = results.get(made_key, 0)
                    print(f"main: made_value is {made_value}")
                    row += "," + str(made_value)
                print(f"main: output row is {row}")
                output_file.write(row + "\n")
            
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Create a matrix from space-delimited files in a directory.")
    parser.add_argument("input_directory", help="Path to the input directory containing space-delimited files.")
    parser.add_argument("output_matrix", help="Path to the output matrix file.")
    args = parser.parse_args()

    main(args.input_directory, args.output_matrix)
