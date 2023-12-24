all:
	snakemake -j 4 -k

html:
	mkdocs build

serve:
	mkdocs serve
