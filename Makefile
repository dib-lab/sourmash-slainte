all:
	snakemake -j 4 -k

sketch:
	snakemake -j 4 -k sketch

html:
	mkdocs build

serve:
	mkdocs serve
