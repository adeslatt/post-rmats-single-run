# post-rmats-single-run
[Kids First Workflow v4](https://github.com/kids-first/kf-rnaseq-workflow)
is an example of running rMATS in the signular that is running it on only one RNA-seq file at a time.


Because the files are treated individually, splice variants, different alternative splice results,
are detected separately and may not be present in every sample, donor, etc.

But to analyze and even classify the separate files and their membership, we have to do the work.


#

1. Cat all the files of the same type together
2. Sort to make a union - this will be the master file
3. Add an ID to the beginning
4. Normalize the files afterwards - which means sorting each of the files
5. and then creating a hash of the master file (the union

To create a normalized file -- 5 separate awk scripts were created.

match_se.awk
match_a3ss.awk
match_a5ss.awk
match_ri.awk
match_mxe.awk

these are each run as follows:

```
awk -f match_se.awk B A
```