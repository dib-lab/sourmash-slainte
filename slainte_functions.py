import glob
import os

__all__ = ['strip_suffix',
           'collect_genomes']

def strip_suffix(x):
    "Remove standard DNA file suffixes to get the filename"
    basename = os.path.basename(x)
    while 1:
        prefix, suffix = os.path.splitext(basename)
        if suffix in ('.fa', '.fna', '.fasta', '.fq', '.fastq', '.gz', '.bz2'):
            basename = prefix
        else:
            break
    return basename

def collect_genomes(genome_path):
    genome_names = {}
    for g in genome_path:
        files = glob.glob(g)
        for filename in files:
            name = strip_suffix(filename)
            assert name not in genome_names, f"duplicate prefix for {filename}"
            genome_names[name] = filename

    return genome_names

