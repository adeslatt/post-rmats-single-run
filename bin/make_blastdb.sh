#!/bin/bash

aafile="*_aa.fa"
blast_end="_blast.db"

echo "aafile is $aafile"

for file in $aafile; do
    name="${file%.SE.coordinates_linear_aa.fa}"
    blastname=$name$blast_end
    title=$name" blast db"
    logfile=$name"_log.txt"

    echo "name is $name"
    echo "blastdb name is $blastname"
    echo "title is $title"
    
    makeblastdb -in $file -dbtype prot -out $blastname
done 
