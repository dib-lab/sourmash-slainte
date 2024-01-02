# this module is used by both Snakefile and mkdocs-macros.
# see 'define_env' for mkdocs-macros hook.
import glob
import os
import csv
from collections import defaultdict
import pandas as pd


__all__ = ['strip_suffix',
           'collect_genomes',
           'load_metagenome_files']

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


def load_metagenome_files(metag_path, samples_csv, debug=False):
    METAGENOME_NAMES=defaultdict(set)
    METAGENOME_FILES=dict()
    with open(samples_csv, 'r', newline='') as sample_fp:
        r = csv.DictReader(sample_fp)

        for row in r:
            # use the 'prefix' column as the prefix for a wildcard
            name = row['name']
            fileglob = os.path.join(metag_path, row['prefix'] + '*')

            files = glob.glob(fileglob)
            pretty_print_files = []
            for f in files:
                pretty_print_files.append(f"'{f}'")
            pretty_print_files = "\n\t" + "\n\t".join(pretty_print_files)
            print(f"for sample '{name}', wildcard '{fileglob}' matches:{pretty_print_files}")
            assert files, fileglob

            METAGENOME_NAMES[name].update(files)

            for filename in files:
                individual_name = strip_suffix(filename)
                assert individual_name not in METAGENOME_FILES, individual_name
                METAGENOME_FILES[individual_name] = filename

    return METAGENOME_NAMES, METAGENOME_FILES


def define_env(env):
    "Hook function for mkdocs"

    @env.macro
    def render_csv(filename):
        if os.path.exists(filename):
            try:
                df = pd.read_csv(filename)
                return df.to_markdown()
            except pd.errors.EmptyDataError:
                pass

        return ""
