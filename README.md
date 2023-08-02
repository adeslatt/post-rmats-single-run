# post-rmats-single-run
[Kids First Workflow v4](https://github.com/kids-first/kf-rnaseq-workflow)
is an example of running [rMATS](https://github.com/Xinglab/rmats-turbo#readme) in the signular that is running it on only one RNA-seq file at a time.

Because the files are treated individually, splice variants, different alternative splice results,
are detected separately and may not be present in every sample, donor, etc.

But to analyze and even classify the separate files and their membership, we have to do the work.

The script `prepareSEfiles.sh` takes the output from supplied single run rMATS analyses and makes `4` matricies:

* `SE.SJC.matrix.txt` - the matrix with the normalized IDs based upon the non-redundant union of all the SE events in supplied files and counts for the skipped exon junctions
* `SE.SJC.w.coordinates.matrix.txt` - containing the coordinates for the exon in question, with upstream and downstream exon coordinates
* `SE.IJC.matrix.txt` - the matrix with the normalized IDS and included junction counts
* `SE.IJC.w.coordinates.matrix.txt` - containing the coordinates for the junction that the counts are concerning.

Also yielding the bed file of all the events that can be loaded as a custom track in the UCSC browser.
* `SE.coordinates.bed`



The coordinates information is the same between the IJC and SJC numbers.  Though somewhat confusing nomenclature from rMATS the IJC counts are most informative in the sense of the positive nature of the counting and confirming information concerning the junction in question.  

The addition of annotation information also helps to keep the user of this information informed regarding the subject of the evidence.

This is very helpful in putting the information together.

## Philosophy

My philosophy - is always an elements of style approach.  

`Input` - `process` - `output`

In the case of preparing the matrices to allow the analysis of alternative splicing with coordinates and counts from multiple samples 

This is done using two awk scripts:

* `match_se.awk`  - this file is used to normalize all the samples put in the matrix so that for each junction their is a uniform ID.  Match uses associated arrays in awk to match the coordinates and if a sample has counts for this junction, they are added, if not, zeros are placed instead creating a non-empty complete matrix - important for analysis.

* [`make_bed_se.awk`](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/make_bed_se.awk) - this file creates from the `SE.coordinates.matrix.txt` an appropriate bed file `SE.coordinates.bed` that may be uploaded as a custom track on the UCSC Genome Browser.

The [`prepareSEfiles.sh`](https://github.com/adeslatt/post-rmats-single-run/blob/main/bin/prepareSEfiles.sh) takes 10 steps to produce the final output.   These steps are documented in the script itself. 

## `TO DO`

* Create similiar matrices for the A3SS, A5SS, RI and MXE data

