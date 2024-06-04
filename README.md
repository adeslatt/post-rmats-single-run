# post-rmats-single-run
[Kids First Workflow v4](https://github.com/kids-first/kf-rnaseq-workflow)
is an example of running [rMATS](https://github.com/Xinglab/rmats-turbo#readme) in the signular that is running it on only one RNA-seq file at a time.

Because the files are treated individually, splice variants, different alternative splice results,
are detected separately and may not be present in every sample, donor, etc.

But to analyze and even classify the separate files and their membership, we have to do the work.

## Software Requirements

All software was installed within a conda environment.   To the greatest extent possible, the installation was done using the packaging offered by Anaconda.  This isn't always possible.  And care must be taken to ensure you are installing the correct version. This is usually the latest version of the tool.  Installing the latest version of the tool isn't always possible because of operating systems.  Just do the best you can :).

* [cpat](https://cpat.readthedocs.io/en/latest/).
  This is the [coding potential assessment tool `(CPAT)`](https://academic.oup.com/nar/article/41/6/e74/2902455) to find open reading frames `(ORFs)`.

  Installed using
  ```bash
  pip3 install CPAT
  ```
 
* [gotranseq](https://github.com/feliixx/gotranseq) - this is used for translating to amino acids the open reading frames found by cpat

* [Bedtools](https://bedtools.readthedocs.io/en/latest/) - this is used to make the bed files for viewing in the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway)

## Inputs

In addition to the output from all of the rMATS runs you want to consider together, you will need the following reference files.

These files are all human relative, because my experiments were done with human samples.  If your samples came from other organisms, you should adjust accordingly.

For cpat and for Bedtools:
* GRCh38.primary_assembly.genome.fa

For cpat alone:
* Human_Hexamer.tsv
* Human_logitModel.RData

## Process

There are five(5) different alternative splicing events that are classified by rMATS:
* SE
* RI
* MXE
* A3SS
* A5SS

The input files and processing files are specific to each of these splicing events and are handled by scripts below individually.

In general, I recommend you make a subdirectory and create symbolic links to your data for each of the separate categories so that they may be analyzed and digested separately.

```bash
mkdir SE_calculate
mkdir RI_calculate
mkdir MXE_calculate
mkdir A3SS_calculate
mkdir A5SS_calculate
```

Then change directory into each of them and make a symbolic link to the files.
For example, assuming that the subdirectories are within the directory where your rMATS output files are 

Please note that the "." means current directory.   This is required for the prepareXXfiles.sh runs (XX== SE, RI, MXE, A3SS, or A5SS)

For the SE files you would do this:
```bash
cd SE_calculate
ln -s ../*SE.MATS.JC.txt .
../../bin/prepareSEfiles.sh .
../../bin/make_individual_SE_files.sh
```

For the RI files, assuming you were in the SE_calculate subdirectory
```bash
cd ../RI_calculate
ln -s ../*RI.MATS.JC.txt .
../../bin/prepareRIfiles.sh .
../../bin/make_individual_RI_files.sh
```

for MXE files, assuming you were in the RI_calculate subdirectory
```bash
cd ../MXE_calculate
ln -s ../*MXE.MATS.JC.txt .
../../bin/prepareMXEfiles.sh .
../../bin/make_individual_MXE_files.sh
```

for A3SS files, assuming you were in the MXE_calculate subdirectory
```bash
cd ../A3SS_calculate
ln -s ../*A3SS.MATS.JC.txt .
../../bin/prepareA3SSfiles.sh .
../../bin/make_individual_A3SS_files.sh
```

for A5SS files, assuming you were in the A3SS_calculate subdirectory
```bash
cd ../A3SS_calculate
ln -s ../*A3SS.MATS.JC.txt .
../../bin/prepareA3SSfiles.sh .
../../bin/make_individual_A5SS_files.sh
```

Now we will look at what is produced by each of these programs.

### Skipped Exon (SE)

The format of the input is as follows:           
```bash                                                                                                                                                                                                                
     col 1 - ID - unique identifier for the skipped exon event 
     col 2 - GeneID - the ENSG identifier                                                                                                                                                                                       
     col 3 - geneSymbol - the text word for the gene
     col 4 - chromsome                                                                                                                                                                                                  
     col 5 - strand
     col 6 - exonStart_0base - the start base coordinate (zero based) for the exon of interest
     col 7 - exonEnd - the end base coordinate for the exon of interest
     col 8 - upstreamES - the start base coordinate (zero based) for the upstream exon
     col 9 - upstreamEE - the end coordinate for the upstream exon
     col 10 - downstreamES - the start base coordinate (zero based) for the downstream exon
     col 11 - downstreamEE - the end coordiante for the downstream exon  
```

The script [prepareSEfiles.sh](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/prepareSEfiles.sh) takes the output from supplied single run rMATS analyses and makes `4` matricies using two awk scripts:
* [match_se.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/match_se.awk)
* [make_bed_se.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/make_bed_se.awk)

* `SE.SJC.matrix.txt` - the matrix with the normalized IDs based upon the non-redundant union of all the SE events in supplied files and counts for the skipped exon junctions
* `SE.SJC.w.coordinates.matrix.txt` - containing the coordinates for the exon in question, with upstream and downstream exon coordinates
* `SE.IJC.matrix.txt` - the matrix with the normalized IDS and included junction counts
* `SE.IJC.w.coordinates.matrix.txt` - containing the coordinates for the junction that the counts are concerning.

Also yielding the bed file of all the events that can be loaded as a custom track in the UCSC browser.
* `SE.coordinates.bed`

### Mutually Excluded Exon (MXE)

The format of the input is as follows:

```bash
     col 1 - ID - unique identifier for the skipped exon event 
     col 2 - GeneID - the ENSG identifier
     col 3 - geneSymbol - the text word for the gene
     col 4 - chromsome
     col 5 - strand  
     col 6 - 1stexonStart_0base - the start base coordinate (zero based) for the exon of interest
     col 7 - 1stexonEnd - the end base coordinate for the exon of interest
     col 8 - 2ndexonStart_0base - start of 2nd exon
     col 9 - 2ndexonEnd
     col 10 - upstreamES - the start base coordinate (zero based) for the upstream exon
     col 11 - upstreamEE - the end coordinate for the upstream exon
     col 12 - downstreamES - the start base coordinate (zero based) for the downstream exon
     col 13 - downstreamEE - the end coordiante for the downstream exon 
```

* If the strand is +, then the inclusion form includes the 1st exon (1stExonStart_0base, 1stExonEnd) and skips the 2nd exon
* If the strand is -, then the inclusion form includes the 2nd exon (2ndExonStart_0base, 2ndExonEnd) and skips the 1st exon

The script [prepareMXEfiles.sh](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/prepareMXEfiles.sh) takes the output from supplied single run rMATS analyses and makes `4` matricies using two awk scripts:

* [match_mxe.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/match_mxe.awk)
* [make_bed_mxe.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/make_bed_mxe.awk)

* `MXE.SJC.matrix.txt` - the matrix with the normalized IDs based upon the non-redundant union of all the SE events in supplied files and counts for the skipped exon junctions
* `MXE.SJC.w.coordinates.matrix.txt` - containing the coordinates for the exon in question, with upstream and downstream exon coordinates
* `MXE.IJC.matrix.txt` - the matrix with the normalized IDS and included junction counts
* `MXE.IJC.w.coordinates.matrix.txt` - containing the coordinates for the junction that the counts are concerning.

Also yielding the bed file of all the events that can be loaded as a custom track in the UCSC browser.
* `MXE.coordinates.bed`

The coordinates information is the same between the IJC and SJC numbers.  Though somewhat confusing nomenclature from rMATS the IJC counts are most informative in the sense of the positive nature of the counting and confirming information concerning the junction in question.  

The addition of annotation information also helps to keep the user of this information informed regarding the subject of the evidence.

This is very helpful in putting the information together.

### Retention Exon (RI)

The format of the input for a Retention Intron

                                                                                                                                                                                                                                 ```bash
     col 1 - ID - unique identifier for the skipped exon event               
     col 2 - GeneID - the ENSG identifier
     col 3 - geneSymbol - the text word for the gene
     col 4 - chromsome
     col 5 - strand
     col 6 - riexonStart_0base - the start base coordinate (zero based) for the exon of interest
     col 7 - riexonEnd - the end base coordinate for the exon of interest 
     col 8 - upstreamES - the start base coordinate (zero based) for the upstream exon
     col 9 - upstreamEE - the end coordinate for the upstream exon
     col 10 - downstreamES - the start base coordinate (zero based) for the downstream exon
     col 11 - downstreamEE - the end coordiante for the downstream exon 
```

The script [prepareRIfiles.sh](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/prepareRIfiles.sh) takes the output from supplied single run rMATS analyses and makes `4` matricies using two awk scripts:
* [match_ri.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/match_ri.awk)
* [make_bed_ri.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/make_bed_ri.awk)

* `RI.SJC.matrix.txt` - the matrix with the normalized IDs based upon the non-redundant union of all the RI events in supplied files and counts for the retention intron junctions
* `RI.SJC.w.coordinates.matrix.txt` - containing the coordinates for the RI exon in question, with upstream and downstream exon coordinates
* `RI.IJC.matrix.txt` - the matrix with the normalized IDS and included junction counts
* `RI.IJC.w.coordinates.matrix.txt` - containing the coordinates for the junction that the counts are concerning.

Also yielding the bed file of all the events that can be loaded as a custom track in the UCSC browser.
* `RI.coordinates.bed`

### Alternative 3' Splice Site (A3SS)

The format of the input is as follows:

```bash
          col 1 - ID 
          col 2 - GeneID 
          col 3 - geneSymbol 
          col 4 - chr 
          col 5 - strand 
          col 6 - longExonStart_0base 
          col 7 - longExonEnd 
          col 8 - shortES 
          col 9 - shortEE
          col 10 - flankingES  
          col 11 - flankingEE
```

The script [prepareA3SSfiles.sh](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/prepareA3SSfiles.sh) takes the output from supplied single run rMATS analyses and makes `4` matricies using two awk scripts:
* [match_a3ss.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/match_a3ss.awk)
* [make_bed_a3ss.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/make_bed_a3ss.awk)

* `A3SS.SJC.matrix.txt` - the matrix with the normalized IDs based upon the non-redundant union of all the A3SS events in supplied files and counts for the retention intron junctions
* `A3SS.SJC.w.coordinates.matrix.txt` - containing the coordinates for the A3SS exon in question, with upstream and downstream exon coordinates
* `A3SS.IJC.matrix.txt` - the matrix with the normalized IDS and included junction counts
* `A3SS.IJC.w.coordinates.matrix.txt` - containing the coordinates for the junction that the counts are concerning.

Also yielding the bed file of all the events that can be loaded as a custom track in the UCSC browser.
* `A3SS.coordinates.bed`

### Alternative 5' Splice Site (A5SS)

The format of the input is as follows:   

```bash
          col 1 - ID  
          col 2 - GeneID 
          col 3 - geneSymbol 
          col 4 - chr  
          col 5 - strand  
          col 6 - longExonStart_0base 
          col 7 - longExonEnd 
          col 8 - shortES 
          col 9 - shortEE
          col 10 - flankingES 
          col 11 - flankingEE
```

The script [prepareA5SSfiles.sh](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/prepareA3SSfiles.sh) takes the output from supplied single run rMATS analyses and makes `4` matricies using two awk scripts:
* [match_a5ss.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/match_a5ss.awk)
* [make_bed_a5ss.awk](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/make_bed_a5ss.awk)

* `A5SS.SJC.matrix.txt` - the matrix with the normalized IDs based upon the non-redundant union of all the A5SS events in supplied files and counts for the retention intron junctions
* `A5SS.SJC.w.coordinates.matrix.txt` - containing the coordinates for the A3SS exon in question, with upstream and downstream exon coordinates
* `A5SS.IJC.matrix.txt` - the matrix with the normalized IDS and included junction counts
* `A5SS.IJC.w.coordinates.matrix.txt` - containing the coordinates for the junction that the counts are concerning.

Also yielding the bed file of all the events that can be loaded as a custom track in the UCSC browser.
* `A5SS.coordinates.bed`

## Getting Domain Matrix Counts

Continuing because our goal is to get domain information and sort through functional differences given our different states for our data collection.  
For this particular study, it is looking for molecular differences between patients who:

* experience TAM and **do not progress** to AML as compared with patients who
* experience TAM and **do progress** to AML

We now get the **protein domains** that are present in these experimental files.

Having run all of the routines for each of the alternative splicing events (SE, RI, MXE, A3SS, A5SS)

Next step to make the search for the exact domain sequence is to linearize these amino acid sequences removing the standard 60 character limitation.
Also run in the same directory.

To perform the grep on each of the files - we need to linearize the AA fasta files.   This is done with the routine `makingAASingleLine.sh`.

Now as a matter of cleanliness and convienence, each of these file types were then put into their own directory - to make the final matrix creation easy.

So within the experiment directory, three subdirectories were made
```bash
mkdir SE_linear
mkdir MXE_linear
mkdir RI_linear
mkdir A3SS_linear
mkdir A5SS_linear
```
Note that you need to specify the directory -- so in this case below you will see `.` which is the current directory.

```bash
 ../bin/makingAASingleLine.sh .
```

And the linear files associated with SE, MXE, and RI moved into their appropriate subdirectories.

## Making protein files

Using [uniprot](https://uniprot.org) we obtain the domain sequences that are part of the Protein.  This way we can arrive at the putative functional differences between those amino acid sequences present in one class of samples versus another.

Manually creating files of interest we create files for our use:

For example, `MYC_human_P01106.txt` looks as follows:

There is no header in the `master` protein sequence file with the domains specified.

Protein_name:Domain_Sequence_Name:Amino-acid Range:Amino-acid Sequence

SO the file for Human MYC then is as follows:

```bash
MYC:9aaTAD:115-123:EMVTELLGG
MYC:Polar_residues:219-249:SPKSCASQDSSAFSPSSDSLLSSTESSPQGS
MYC:Disordered:219-310:SPKSCASQDSSAFSPSSDSLLSSTESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEEEIDVVSVEKRQAPGKRSESGSPSAGGHSKPPHSPL
MYC:bHLH:369-421:VKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILSV
MYC:Leucine_zipper:428-449:LISEEDLLRKRREQLKHKLEQL
```

By default I name the file with the Uniprot Identifier -- just as a convention to help my future self.


## Counting the Domain hits

Now in another directory here named `protein_domain_counts` we create a subdirectory for each of our proteins and splicing types:
With MYC the following sub directories are made:
```bash
mkdir protein_domain_counts/MYC_MXE
mkdir protein_domain_counts/MYC_SE
mkdir protein_domain_counts/MYC_RI
mkdir protein_domain_counts/MYC_A3SS
mkdir protein_domain_counts/MYC_A5SS
```

We change directory into each of them and count, for example:
```bash
cd protein_domain_counts/MYC_SE
awk -v experiment_dir="/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/paired.TAM.AMLv2/SE_linear/" -f /Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/bin/process_domains.awk "/Users/annedeslattesmays/Desktop/projects/post-rmats-single-run/protein_aa/MYC_human_P01106.txt"
```

Now we have each of the experiments and the read counts - but to analyze better a matrix would serve the rule.
Make a directory and create the appropriately named matrix

```bash
cd ../..
mkdir protein_domain_matrices
cd protein_domain_matrices
```
```bash
python ../bin/make_protein_counts_matrix.py ../protein_domain_counts/MYC_SE/ MYC_SE_protein_domain_counts.csv > debug.txt
```

```bash
Protein,Domain_Name,AA_position,Domain_Sequence,PAUVKY-03A-01R.SE.coordinates_linear_aa.fa_results.txt,PAUVKY-40A-01R.SE.coordinates_linear_aa.fa_results.txt,PAWHSD-03A-01R.SE.coordinates_linear_aa.fa_results.txt,PAWHSD-40A-01R.SE.coordinates_linear_aa.fa_results.txt,PAWSNZ-03A-01R.SE.coordinates_linear_aa.fa_results.txt,PAWSNZ-40A-01R.SE.coordinates_linear_aa.fa_results.txt,_1_PAUTLA-03A-01R.SE.coordinates_linear_aa.fa_results.txt,_1_PAUTLA-40A-01R.SE.coordinates_linear_aa.fa_results.txt,_1_PAVUDU-03A-01R.SE.coordinates_linear_aa.fa_results.txt,_1_PAVUDU-40A-01R.SE.coordinates_linear_aa.fa_results.txt
MYC,9aaTAD,115-123,EMVTELLGG,0,0,2,0,1,0,1,0,2,0
MYC,Polar_residues,219-249,SPKSCASQDSSAFSPSSDSLLSSTESSPQGS,0,0,2,0,0,0,0,0,2,0
MYC,Disordered,219-310,SPKSCASQDSSAFSPSSDSLLSSTESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEEEIDVVSVEKRQAPGKRSESGSPSAGGHSKPPHSPL,0,0,0,0,0,0,0,0,0,0
MYC,bHLH,369-421,VKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILSV,0,0,0,0,0,0,0,0,0,0
MYC,Leucine_zipper,428-449,LISEEDLLRKRREQLKHKLEQL,0,0,0,0,0,0,0,0,0,0
```

Now we can analyze.
