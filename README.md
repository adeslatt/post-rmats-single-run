# post-rmats-single-run
[Kids First Workflow v4](https://github.com/kids-first/kf-rnaseq-workflow)
is an example of running [rMATS](https://github.com/Xinglab/rmats-turbo#readme) in the signular that is running it on only one RNA-seq file at a time.

Because the files are treated individually, splice variants, different alternative splice results,
are detected separately and may not be present in every sample, donor, etc.

But to analyze and even classify the separate files and their membership, we have to do the work.

## Inputs

Output from an rMATS run

## Process

### Skipped Exon (SE)

From rMATS documentation, we learn the SE event: exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE
The inclusion form (IJC) includes the target exon (exonStart_0base, exonEnd)

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

From the rMATS documentation we learn that the MXE event: MXE: 1stExonStart_0base 1stExonEnd 2ndExonStart_0base 2ndExonEnd upstreamES upstreamEE downstreamES downstreamEE
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

## Getting Domain Matrix Counts

Continuing because our goal is to get domain information and sort through functional differences given our different states for our data collection.  
For this particular study, it is looking for molecular differences between patients who:

* experience TAM and **do not progress** to AML as compared with patients who
* experience TAM and **do progress** to AML

We now get the **protein domains** that are present in these experimental files.

Using the coordinate files we obtain through a series of steps:

* make bed files (via 'make_se_bed.awk', 'make_ri_bed.awk', 'make_mxe_bed.awk')
* obtain the sequence data covered by the regions outlined by the measured results (via 'bedtools getfasta')
* call the protein open reading frames covered by these sequence data (via 'cpat.py')
* translate these protein ORFs to amino acid sequence (via 'gotranseq')

This is obtained with the execution of `make_individual_files.sh'

Run from the experimental data directory:
```bash
../../bin/make_individual_files.sh
```

Next step to make the search for the exact domain sequence is to linearize these amino acid sequences removing the standard 60 character limitation.
Also run in the same directory.

```bash
 ../bin/makingAASingleLine.sh .
```

Now as a matter of cleanliness and convienence, each of these file types were then put into their own directory - to make the final matrix creation easy.

So within the experiment directory, three subdirectories were made
```bash
mkdir SE_linear
mkdir MXE_linear
mkdir RI_linear
```
And the linear files associated with SE, MXE, and RI moved into their appropriate subdirectories.

## Making protein files

Using [uniprot](https://uniprot.org) we obtain the domain sequences that are part of the Protein.  This way we can arrive at the putative functional differences between those amino acid sequences present in one class of samples versus another.

Manually creating files of interest we create files for our use:

For example, `MYC_human_P01106.txt` looks as follows:
```bash
MYC:9aaTAD:115-123:EMVTELLGG
MYC:Polar_residues:219-249:SPKSCASQDSSAFSPSSDSLLSSTESSPQGS
MYC:Disordered:219-310:SPKSCASQDSSAFSPSSDSLLSSTESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEEEIDVVSVEKRQAPGKRSESGSPSAGGHSKPPHSPL
MYC:bHLH:369-421:VKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILSV
MYC:Leucine_zipper:428-449:LISEEDLLRKRREQLKHKLEQL
```
The file is a table separated by ':' (colon).  The fields are in order:
* Protein name
* Domain description
* Amino acid positions
* Amino acid sequence

## Counting the Domain hits

Now in another directory here named `protein_domain_counts` we create a subdirectory for each of our proteins and splicing types:
With MYC the following sub directories are made:
```bash
mkdir protein_domain_counts/MYC_MXE
mkdir protein_domain_counts/MYC_SE
mkdir protein_domain_counts/MYC_RI
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

## Philosophy

My philosophy - is always an elements of style approach.  

`Input` - `process` - `output`

In the case of preparing the matrices to allow the analysis of alternative splicing with coordinates and counts from multiple samples 

This is done using two awk scripts:

* [`match_se.awk`](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/match_se.awk)  - this file is used to normalize all the samples put in the matrix so that for each junction their is a uniform ID.  Match uses associated arrays in awk to match the coordinates and if a sample has counts for this junction, they are added, if not, zeros are placed instead creating a non-empty complete matrix - important for analysis.

* [`make_bed_se.awk`](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/make_bed_se.awk) - this file creates from the `SE.coordinates.matrix.txt` an appropriate bed file `SE.coordinates.bed` that may be uploaded as a custom track on the UCSC Genome Browser.

The [`prepareSEfiles.sh`](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/prepareSEfiles.sh) takes 10 steps to produce the final output.   These steps are documented in the script itself. 

## `TO DO`

* Create similiar matrices for the A3SS, A5SS
* Turn into a proper Nextflow workflow to run on a platform - Lifebit's Data Science Suite Module (previously CloudOS).

