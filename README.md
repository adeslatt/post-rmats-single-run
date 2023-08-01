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

The coordinates information is the same between the IJC and SJC numbers.  Though somewhat confusing nomenclature from rMATS the IJC counts are most informative in the sense of the positive nature of the counting and confirming information concerning the junction in question.  

The addition of annotation information also helps to keep the user of this information informed regarding the subject of the evidence.

This is very helpful in putting the information together.

## `TO DO`

* Create similiar matrices for the A3SS, A5SS, RI and MXE data

