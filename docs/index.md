# Welcome to sourmash slainte

slainte does basic comparisons and analyses of metagenomes and genomes
using [sourmash](https://sourmash.readthedocs.io/).

[Metagenome comparisons](metag_compare.md) - which metagenomes are similar?

{% if display_genomes %}
[Genome comparisons](genome_compare.md) - which genomes are similar?

[Metagenomes against genomes](metag_x_genomes.md) - which genomes are in which metagenomes?
{% endif %}

{% if run_gather %}
[Metagenome taxonomy summary](metag_tax.md) - a taxonomic breakdown across all metagenomes.
{% endif %}

---

<!-- [Config / macros information](macros_info.md) -->

---

{% if run_gather %}

## Databases and taxonomy

Search databases used for `sourmash gather`:
{% for db in databases %}
* {{ db }}
{% endfor %}

Taxonomy used for `sourmash tax metagenome`:
{% for tax in taxonomies %}
* {{ tax }}
{% endfor %}

{% else %}

No metagenome gather or taxonomic analysis is being performed;
`run_gather` setting in `config.yml` is False.

{% endif %}

## Metagenome sample information

See the [sample x datafile](metag_sample_check.md) comparison to
verify which samples contain which data.
