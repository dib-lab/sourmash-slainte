# sourmash-slainte

## Quickstart

```
srun -p high2 --time=72:00:00 --nodes=1 --cpus-per-task 128 --mem 120GB --pty /bin/bash
```

```
mamba env create -n slainte -f environment.yml 
```

```
snakemake --cores 128 --resources gather=2 -k
```

estimate ~20 GB per gather?
fastgather doesn't use more memory with threads

