# sourmash-slainte

sourmash-slainte is a snakemake workflow that takes in a collection of
metagenomes and (optionally) a set of genome sequences (e.g. reference
genomes, isolates, SAGs, or MAGs).  It then performs a number of
default comparisons. sourmash-slainte is intended to be a good
starting template for projects that need to work with many large
metagenomes.

The current list of default analyses (Jan 2024) is:
* all-by-all comparison of metagenomes
* all-by-all comparison of genomes
* genome detection within metagenomes
* genomic composition of metagenomes & associated taxonomy

## Quickstart

Allocate computational resources, here 128 CPUs and 120 GB of RAM for
3 days via Slurm:

```
srun -p high2 --time=72:00:00 --nodes=1 --cpus-per-task 128 --mem 120GB --pty /bin/bash
```

Create a conda environment with the necessary software:

```
mamba env create -n slainte -f environment.yml 
```

Run!
```
snakemake --cores 128 --resources gather=2 -k
```

CTB to add:
* some minimal configuration!
* estimate ~20 GB per gather?
* n.b. fastgather doesn't use more memory with threads
