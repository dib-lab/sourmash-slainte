# Genome comparisons

{% if display_genomes %}

## ANI comparison between genomes

These plots are the output of `sourmash compare --containment --ani`
[(link)](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-compare-compare-many-signatures)
and `sourmash plot`
[(link)](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-plot-cluster-and-visualize-comparisons-of-many-signatures).

The comparison metric displayed here is a k-mer-based estimation of
[the Average Nucleotide Identity](https://en.wikipedia.org/wiki/Bacterial_genome#Genome_comparisons). Genome
pairs with similar nucleotide content will cluster together.

### k=21

[![](outputs/genome_compare.21.ani.matrix.png)](outputs/genome_compare.21.ani.matrix.png)

### k=31

[![](outputs/genome_compare.31.ani.matrix.png)](outputs/genome_compare.31.ani.matrix.png)

### k=51

[![](outputs/genome_compare.51.ani.matrix.png)](outputs/genome_compare.51.ani.matrix.png)

{% else %}

(No individual query genomes provided, so there's nothing to show.)

{% endif %}
