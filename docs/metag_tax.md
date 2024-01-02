# Metagenome taxonomic summary

These tables are created by running `sourmash gather` followed by
`sourmash tax metagenome` on each metagenome sample.

Search databases used for `sourmash gather`:
{% for db in databases %}
* {{ db }}
{% endfor %}

Taxonomy used for `sourmash tax metagenome`:
{% for tax in taxonomies %}
* {{ tax }}
{% endfor %}

## Taxonomic summary

{{ read_csv('outputs/metag_gather/metag.21.kreport.csv') }}

[(Download CSV file)](outputs/metag_gather/metag.21.kreport.csv)
