#!/bin/bash

allA3SSoriginalFiles="*.filtered.A3SS.MATS.JC.txt"
rmats="rMATS_"
A3SS="A3SS"
description="A3SS Events for DS-AML"
track_name="track name="
description=" description="
coordinates_txt_end=".coordinates.txt"
coordinates_bed_end=".coordinates.bed"
coordinates_fasta_end=".coordinates.fasta"
cpat_output=".cpat_output"
cpat_error=".cpat_error"
orf_seqs_end=".ORF_seqs.fa"
orf_linear_aa_end="_linear_aa.fa"
orf_seqs_aa_end=".ORF_seqs_aa.fa"

underscore="_"
space=" "
period="."

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "SCRIPT_DIR $SCRIPT_DIR"

# Construct the relative path to the gotranseq directory
TRANSEQ_DIR="$SCRIPT_DIR/../../gotranseq"
echo "TRANSEQ_DIR is $TRANSEQ_DIR"

# make the individual sample bed files
for file in $allA3SSoriginalFiles; do
    name="${file%_RBS_withJunctionsOnGenome_dupsFlagged_r1.filtered.A3SS.MATS.JC.txt}"
    header=$track_name$rmats$name$underscore$A3SS$description$rmats$name$underscore$A3SS

    echo "name is $name"
    echo "header is $header"

    txt_filename=$name$period$A3SS$coordinates_txt_end
    bed_filename=$name$period$A3SS$coordinates_bed_end
    fasta_filename=$name$period$A3SS$coordinates_fasta_end
    orf_filename="${txt_filename%.txt}"
    cpat_output_filename=$name$period$A3SS$cpat_output
    cpat_error_filename=$name$period$A3SS$cpat_error
    orf_seqs_filename=$orf_filename$orf_seqs_end
    orf_seqs_aa_filename=$orf_filename$orf_seqs_aa_end
    orf_linear_aa_filename=$orf_filename$orf_linear_aa_end
    
    echo "txt_filename is $orf_filename"
    echo "bed_filename is $bed_filename"
    echo "fasta_filename is $fasta_filename"
    echo "cpat_output_filename is $cpat_output_filename"
    echo "cpat_error_filename is $cpat_error_filename"
    echo "orf_seqs_aa_filename is $orf_seqs_aa_filename"
    echo "orf_linear_aa_filename is $orf_linear_aa_filename"
    
    cut -f 1-11 $file > $txt_filename

    echo $header > $bed_filename
    
    # old absolute path
    #     awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_a3ss.awk $file >> $bed_filename
    # new relative path
    gawk -f "$SCRIPT_DIR/make_bed_a3ss.awk" $file >> $bed_filename
    
    bedtools getfasta -rna -fi GRCh38.primary_assembly.genome.fa -bed $bed_filename > $fasta_filename
    cpat -x $human_hexamer -d $human_logitmodel -x Human_Hexamer.tsv -d Human_logitModel.RData.gz -g $fasta_filename --min-orf=50 --top-orf=50 -o $orf_filename 1> $cpat_output_filename 2> $cpat_error_filename
    "$TRANSEQ_DIR/gotranseq" --sequence $orf_seqs_filename -o $orf_seqs_aa_filename -f 1
    
done
