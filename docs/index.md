# Welcome to sourmash slainte

slainte does basic comparisons and analyses of metagenomes and genomes
using [sourmash](https://sourmash.readthedocs.io/).

[Metagenome comparisons](metag_compare.md) - which metagenomes are similar?

{% if display_genomes %}
[Genome comparisons](genome_compare.md) - which genomes are similar?

[Metagenomes against genomes](metag_x_genomes.md) - which genomes are in which metagenomes?
{% endif %}

[Metagenome taxonomy summary](metag_tax.md) - a taxonomic breakdown across all metagenomes.

---

<!-- [Config / macros information](macros_info.md) -->

---

## Databases and taxonomy

Search databases used for `sourmash gather`:
{% for db in databases %}
* {{ db }}
{% endfor %}

Taxonomy used for `sourmash tax metagenome`:
{% for tax in taxonomies %}
* {{ tax }}
{% endfor %}

## Sample information

See the [sample x datafile](metag_sample_check.md) comparison for
which samples contain which data.
