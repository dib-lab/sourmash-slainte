all:
	snakemake -j 4 -k

sketch:
	snakemake -j 4 -k sketch

cleanall:
	rm -fr sketches outputs

html:
	mkdocs build

serve:
	mkdocs serve --watch Snakefile --watch config.yml
