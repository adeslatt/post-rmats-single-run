#!/usr/bin/env nextflow

params.allSEoriginalFiles = "*.filtered.SE.MATS.JC.txt"
params.allMXEoriginalFiles = "*.filtered.MXE.MATS.JC.txt"
params.allRIoriginalFiles = "*.filtered.RI.MATS.JC.txt"
params.rmats = "rMATS_"
params.RI = "RI"
params.SE = "SE"
params.MXE = "MXE"
params.description = " Events for DS-AML"
params.track_name = "track name="
params.description = " description="
params.coordinates_txt_end = ".coordinates.txt"
params.coordinates_bed_end = ".coordinates.bed"
params.coordinates_fasta_end = ".coordinates.fasta"
params.cpat_output = ".cpat_output"
params.cpat_error = ".cpat_error"
params.orf_seqs_end = ".ORF_seqs.fa"
params.orf_linear_aa_end = "_linear_aa.fa"
params.orf_seqs_aa_end = ".ORF_seqs_aa.fa"

// Define channels for input files
inputSE = file(params.allSEoriginalFiles)
inputMXE = file(params.allMXEoriginalFiles)
inputRI = file(params.allRIoriginalFiles)

// Define the process to work on SE files
process processSE {
    input:
    file(fileSE) from inputSE

    output:
    file('*.txt') into txt_files
    file('*.bed') into bed_files
    file('*.fasta') into fasta_files
    file('*.cpat_output') into cpat_output_files
    file('*.cpat_error') into cpat_error_files

    script:
    """
    name = fileSE.baseName.replaceAll('_RBS_withJunctionsOnGenome_dupsFlagged_r1.filtered.SE.MATS.JC', '')
    header = '${params.track_name}${params.rmats}${name}_${params.SE}${params.description}${params.rmats}${name}_${params.SE}'

    echo "name is $name"
    echo "header is $header"

    txt_filename = "${name}.${params.SE}${params.coordinates_txt_end}"
    bed_filename = "${name}.${params.SE}${params.coordinates_bed_end}"
    fasta_filename = "${name}.${params.SE}${params.coordinates_fasta_end}"
    orf_filename = "${name}${params.SE}${params.orf_seqs_end}"
    cpat_output_filename = "${name}.${params.SE}${params.cpat_output}"
    cpat_error_filename = "${name}.${params.SE}${params.cpat_error}"
    orf_seqs_aa_filename = "${name}${params.SE}${params.orf_seqs_aa_end}"
    orf_linear_aa_filename = "${name}${params.SE}${params.orf_linear_aa_end}"
    
    echo "txt_filename is $orf_filename"
    echo "bed_filename is $bed_filename"
    echo "fasta_filename is $fasta_filename"
    echo "cpat_output_filename is $cpat_output_filename"
    echo "cpat_error_filename is $cpat_error_filename"
    echo "orf_seqs_aa_filename is $orf_seqs_aa_filename"
    echo "orf_linear_aa_filename is $orf_linear_aa_filename"
    
    cut -f 1-11 $fileSE > $txt_filename

    echo $header > $bed_filename
    
    awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_se.awk $txt_filename >> $bed_filename
    bedtools getfasta -rna -fi ../../singleCellLongReadAnalysis/data/BC_ranked_isoforms/GRCh38.primary_assembly.genome.fa -bed $bed_filename > $fasta_filename
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/cpat:addr cpat.py -x $human_hexamer -d $human_logitmodel -x Human_Hexamer.tsv -d Human_logitModel.RData.gz -g $fasta_filename --min-orf=50 --top-orf=50 -o $orf_filename 1> $cpat_output_filename 2> $cpat_error_filename
    ../../gotranseq/gotranseq --sequence $orf_seqs_filename -o $orf_seqs_aa_filename -f 1
    """
}

// Define the process to work on RI files
process processRI {
    input:
    file(fileRI) from inputRI

    output:
    file('*.txt') into txt_files
    file('*.bed') into bed_files
    file('*.fasta') into fasta_files
    file('*.cpat_output') into cpat_output_files
    file('*.cpat_error') into cpat_error_files

    script:
    """
    name = fileRI.baseName.replaceAll('_RBS_withJunctionsOnGenome_dupsFlagged_r1.filtered.RI.MATS.JC', '')
    header = '${params.track_name}${params.rmats}${name}_${params.RI}${params.description}${params.rmats}${name}_${params.RI}'

    echo "name is $name"
    echo "header is $header"

    txt_filename = "${name}.${params.RI}${params.coordinates_txt_end}"
    bed_filename = "${name}.${params.RI}${params.coordinates_bed_end}"
    fasta_filename = "${name}.${params.RI}${params.coordinates_fasta_end}"
    orf_filename = "${name}${params.RI}${params.orf_seqs_end}"
    cpat_output_filename = "${name}.${params.RI}${params.cpat_output}"
    cpat_error_filename = "${name}.${params.RI}${params.cpat_error}"
    orf_seqs_aa_filename = "${name}${params.RI}${params.orf_seqs_aa_end}"
    orf_linear_aa_filename = "${name}${params.RI}${params.orf_linear_aa_end}"
    
    echo "txt_filename is $orf_filename"
    echo "bed_filename is $bed_filename"
    echo "fasta_filename is $fasta_filename"
    echo "cpat_output_filename is $cpat_output_filename"
    echo "cpat_error_filename is $cpat_error_filename"
    echo "orf_seqs_aa_filename is $orf_seqs_aa_filename"
    echo "orf_linear_aa_filename is $orf_linear_aa_filename"
    
    cut -f 1-11 $fileRI > $txt_filename

    echo $header > $bed_filename
    
    awk -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/make_bed_ri.awk $txt_filename >> $bed_filename
    bedtools getfasta -rna -fi ../../singleCellLongReadAnalysis/data/BC_ranked_isoforms/GRCh38.primary_assembly.genome.fa -bed $bed_filename > $fasta_filename
    docker run -v $PWD:$PWD -w $PWD -it gsheynkmanlab/cpat:addr cpat.py -x $human_hexamer -d $human_logitmodel -x Human_Hexamer.tsv -d Human_logitModel.RData.gz -g $fasta_filename --min-orf=50 --top-orf=50 -o $orf_filename 1> $cpat_output_filename 2> $cpat_error_filename
    
    ../../gotranseq/gotranseq --sequence $orf_seqs_filename -o $
