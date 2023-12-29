# Metagenome comparisons

## Angular similarity (using abundance information)

These plots are the output of `sourmash compare`
[(link)](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-compare-compare-many-signatures)
and `sourmash plot`
[(link)](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-plot-cluster-and-visualize-comparisons-of-many-signatures).

The comparison metric displayed here is 1 - the
[angular distance](https://en.wikipedia.org/wiki/Cosine_similarity#Angular_distance_and_similarity). Metagenome pairs with similar abundance-weighted nucleotide content will cluster together.

### k=21

[![](outputs/metag_compare.21.abund.matrix.png)](outputs/metag_compare.21.abund.matrix.png)

### k=31

[![](outputs/metag_compare.31.abund.matrix.png)](outputs/metag_compare.31.abund.matrix.png)

### k=51

[![](outputs/metag_compare.51.abund.matrix.png)](outputs/metag_compare.51.abund.matrix.png)

## Jaccard similarity (flat - no abundance used)

These plots are the output of `sourmash compare --ignore-abundance`
[(link)](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-compare-compare-many-signatures)
and `sourmash plot`
[(link)](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-plot-cluster-and-visualize-comparisons-of-many-signatures).

The comparison metric displayed here is the
[Jaccard index or Jaccard similarity](https://en.wikipedia.org/wiki/Jaccard_index). Metagenome
pairs with similar (unweighted) nucleotide content will cluster
together. These comparisons use presence/absence of genomiccontent;
unlike the angular similarity plots above, increased or decreased
abundance of genomic content does not change the relationships.

### k=21

[![](outputs/metag_compare.21.flat.matrix.png)](outputs/metag_compare.21.flat.matrix.png)

### k=31

[![](outputs/metag_compare.31.flat.matrix.png)](outputs/metag_compare.31.flat.matrix.png)

### k=51

[![](outputs/metag_compare.51.flat.matrix.png)](outputs/metag_compare.51.flat.matrix.png)

