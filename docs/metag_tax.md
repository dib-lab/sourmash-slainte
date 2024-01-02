# Metagenome taxonomic summary

{% if run_gather %}

These tables are created by running `sourmash gather` followed by
`sourmash tax metagenome` on each metagenome sample.

## Taxonomic summary

{{ render_csv('outputs/metag_gather/metag.21.kreport.csv') }}

[(Download CSV file)](outputs/metag_gather/metag.21.kreport.csv)

{% else %}

No taxonomic analysis run - `run_gather` is set to False in the config file.

{% endif %}
