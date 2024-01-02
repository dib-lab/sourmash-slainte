# Metagenomes vs genomes

{% if display_genomes %}

## Detection of genomes in metagenomes

The plots below show "detection", or k-mer containment, of provided
reference genomes in the provided metagenomes - you can read more
details
[in the sourmash documentation](https://sourmash.readthedocs.io/en/latest/faq.html#mapping-reads-to-reference-vs-k-mer-detection-or-containment). Intuitively,
a high detection means that reads from that metagenome will cover most
of the genome.

High detection (lighter colors) indicates that genome is present in
that metagenome; detection between 0% and ~75% (red) suggests that a
genome from the same species is present in the metagenome; and no
detection (black) suggests that species is not present.

The plots are currently constructed by running
[manysearch](https://github.com/sourmash-bio/sourmash_plugin_branchwater/tree/main/doc#running-manysearch)
from the
[branchwater plugin for sourmash](https://github.com/sourmash-bio/sourmash_plugin_branchwater).

### k=21

[![](outputs/metag.x.genomes.21.manysearch.png)](outputs/metag.x.genomes.21.manysearch.png)

Files:

* [manysearch output](outputs/metag.x.genomes.21.manysearch.raw.csv)
* [summary matrix output](outputs/metag.x.genomes.21.manysearch.summary.csv)

### k=31

[![](outputs/metag.x.genomes.31.manysearch.png)](outputs/metag.x.genomes.31.manysearch.png)

Files:

* [manysearch output](outputs/metag.x.genomes.31.manysearch.raw.csv)
* [summary matrix output](outputs/metag.x.genomes.31.manysearch.summary.csv)

### k=51

[![](outputs/metag.x.genomes.51.manysearch.png)](outputs/metag.x.genomes.51.manysearch.png)

Files:

* [manysearch output](outputs/metag.x.genomes.51.manysearch.raw.csv)
* [summary matrix output](outputs/metag.x.genomes.51.manysearch.summary.csv)

{% else %}

(No individual query genomes provided, so there's nothing to show.)

{% endif %}
