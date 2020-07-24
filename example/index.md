

# metatree

## Example

### Input:
The test files can be found under the [example/input](https://github.com/aaronmussig/metatree/example/input) directory.

```shell script
metatree batchfile.tsv /tmp/metatree_out taxonomy.tsv p__Altiarchaeota 5
```

### Output:

##### Console:
```
[2020-07-24 11:01:36] INFO: metatree v0.0.1
[2020-07-24 11:01:36] INFO: metatree batchfile.tsv /tmp/output taxonomy.tsv p__Altiarchaeota 5
[2020-07-24 11:01:36] INFO: Rooting trees using GenomeTreeTk v0.1.6
100%|██████████████████████████████| 5/5 [00:01<00:00,  4.59it/s]
[2020-07-24 11:01:37] INFO: Decorating trees using Phylorank v0.1.3
100%|██████████████████████████████| 5/5 [00:07<00:00,  1.45s/it]
[2020-07-24 11:01:45] INFO: Robinson-Foulds metrics will only consider those 1,178 taxa which are common between ALL trees.
[2020-07-24 11:01:45] INFO: Calculating Robinson-Foulds distances.
100%|██████████████████████████████| 10/10 [00:03<00:00,  3.32it/s]
[2020-07-24 11:01:48] INFO: Pairwise Robinson-Foulds distances written to: /tmp/output/results/robinson_foulds_common_taxa/rf_common_taxa.tsv
[2020-07-24 11:01:48] INFO: Calculating Robinson-Foulds distances.
100%|██████████████████████████████| 10/10 [00:00<00:00, 10.49it/s]
[2020-07-24 11:01:49] INFO: Pairwise Robinson-Foulds distances written to: /tmp/output/results/robinson_foulds_all_taxa/rf_all_taxa.tsv
[2020-07-24 11:01:49] INFO: Writing pairwise Robinson-Foulds distances for common taxa to: /tmp/output/results/robinson_foulds_common_taxa
[2020-07-24 11:01:50] INFO: Writing pairwise Robinson-Foulds distances for all taxa to: /tmp/output/results/robinson_foulds_all_taxa
[2020-07-24 11:01:52] INFO: Done.
```

##### Files:
Note that this is not a complete list of all output files, only those which are
of primary interest.

**Robinson-Foulds distances between trees:**

<img src="https://github.com/aaronmussig/metatree/raw/master/example/output/rf_normed_heatmap.svg" alt="Normalised RF heatmap" width="500" height="500">

**PhyloRank polyphyly overview:**

![PhyloRank polyphyly overview](https://github.com/aaronmussig/metatree/raw/master/example/output/tree_comparison_legend.svg "PhyloRank polyphyly overview")
