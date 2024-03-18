#!/bin/bash

allSEoriginalFiles="*.filtered.SE.MATS.JC.txt"
allMXEoriginalFiles="*.filtered.MXE.MATS.JC.txt"
allRIoriginalFiles="*.filtered.RI.MATS.JC.txt"
rmats="rMATS_"
RI="RI"
SE="SE"
MXE="MXE"
description=" Events for DS-AML"
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

# make the individual sample bed files
for file in $allSEoriginalFiles; do
    name="${file%_RBS_withJunctionsOnGenome_dupsFlagged_r1.filtered.SE.MATS.JC.txt}"
    header=$track_name$rmats$name$underscore$SE$description$rmats$name$underscore$SE

    echo "name is $name"
    echo "header is $header"

    txt_filename=$name$period$SE$coordinates_txt_end
    bed_filename=$name$period$SE$coordinates_bed_end
    fasta_filename=$name$period$SE$coordinates_fasta_end
    orf_filename="${txt_filename%.txt}"
    cpat_output_filename=$name$period$SE$cpat_output
    cpat_error_filename=$name$period$SE$cpat_error
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
    
    awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_se.awk $txt_filename >> $bed_filename
    bedtools getfasta -rna -fi ../../singleCellLongReadAnalysis/data/BC_ranked_isoforms/GRCh38.primary_assembly.genome.fa -bed $bed_filename > $fasta_filename
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/cpat:addr cpat.py -x $human_hexamer -d $human_logitmodel -x Human_Hexamer.tsv -d Human_logitModel.RData.gz -g $fasta_filename --min-orf=50 --top-orf=50 -o $orf_filename 1> $cpat_output_filename 2> $cpat_error_filename
    ../../gotranseq/gotranseq --sequence $orf_seqs_filename -o $orf_seqs_aa_filename -f 1

    awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' "$file" > $orf_linear_aa_filename

done

for file in $allRIoriginalFiles; do
    name="${file%_RBS_withJunctionsOnGenome_dupsFlagged_r1.filtered.RI.MATS.JC.txt}"
    header=$track_name$rmats$name$underscore$RI$description$rmats$name$underscore$RI

    echo "name is $name"
    echo "header is $header"

    txt_filename=$name$period$RI$coordinates_txt_end
    bed_filename=$name$period$RI$coordinates_bed_end
    fasta_filename=$name$period$RI$coordinates_fasta_end
    orf_filename="${txt_filename%.text}"
    cpat_output_filename=$name$period$RI$cpat_output
    cpat_error_filename=$name$period$RI$cpat_error
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
    
    awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_ri.awk $txt_filename >> $bed_filename
    bedtools getfasta -rna -fi ../../singleCellLongReadAnalysis/data/BC_ranked_isoforms/GRCh38.primary_assembly.genome.fa -bed $bed_filename > $fasta_filename
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/cpat:addr cpat.py -x $human_hexamer -d $human_logitmodel -x Human_Hexamer.tsv -d Human_logitModel.RData.gz -g $fasta_filename --min-orf=50 --top-orf=50 -o $orf_filename 1> $cpat_output_filename 2> $cpat_error_filename
    
    ../../gotranseq/gotranseq --sequence $orf_seqs_filename -o $orf_seqs_aa_filename -f 1

    awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' "$file" > $orf_linear_aa_filename

done

for file in $allMXEoriginalFiles; do
    name="${file%_RBS_withJunctionsOnGenome_dupsFlagged_r1.filtered.MXE.MATS.JC.txt}"
    header=$track_name$rmats$name$underscore$MXE$description$rmats$name$underscore$MXE

    echo "name is $name"
    echo "header is $header"

    txt_filename=$name$period$MXE$coordinates_txt_end
    bed_filename=$name$period$MXE$coordinates_bed_end
    txt_filename=$name$period$MXE$coordinates_txt_end
    bed_filename=$name$period$MXE$coordinates_bed_end
    fasta_filename=$name$period$MXE$coordinates_fasta_end
    orf_filename="${txt_filename%.text}"
    cpat_output_filename=$name$period$MXE$cpat_output
    cpat_error_filename=$name$period$MXE$cpat_error    
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
    
    cut -f 1-13 $file > $txt_filename

    echo $header > $bed_filename
    
    awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_mxe.awk $txt_filename >> $bed_filename
    bedtools getfasta -rna -fi ../../singleCellLongReadAnalysis/data/BC_ranked_isoforms/GRCh38.primary_assembly.genome.fa -bed $bed_filename > $fasta_filename
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/cpat:addr cpat.py -x $human_hexamer -d $human_logitmodel -x Human_Hexamer.tsv -d Human_logitModel.RData.gz -g $fasta_filename --min-orf=50 --top-orf=50 -o $orf_filename 1> $cpat_output_filename 2> $cpat_error_filename
    ../../gotranseq/gotranseq --sequence $orf_seqs_filename -o $orf_seqs_aa_filename -f 1

    awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' "$file" > $orf_linear_aa_filename

done


    

