import glob, os

configfile: "config.yml"

wildcard_constraints:
    name='[^./]+'

def strip_suffix(x):
    basename = os.path.basename(x)
    while 1:
        prefix, suffix = os.path.splitext(basename)
        if suffix in ('.fa', '.fasta', '.fq', '.fastq', '.gz', '.bz2'):
            basename = prefix
        else:
            break
    print(x, basename)
    return basename

METAGENOME_FILE_PATTERN = config['metagenomes']
GENOME_FILE_PATTERN = config['genomes']

GENOME_NAMES = {}
for g in GENOME_FILE_PATTERN:
    files = glob.glob(g)
    for filename in files:
        name = strip_suffix(filename)
        assert name not in GENOME_NAMES, f"duplicate prefix for {filename}"
        GENOME_NAMES[name] = filename

print(f"Found {len(GENOME_NAMES)} genome files.")

METAGENOME_NAMES={}
for g in METAGENOME_FILE_PATTERN:
    files = glob.glob(g)
    for filename in files:
        name = strip_suffix(filename)
        assert name not in METAGENOME_NAMES, f"duplicate prefix for {filename}"
        METAGENOME_NAMES[name] = filename

print(f"Found {len(METAGENOME_NAMES)} metagenome files.")


rule all:
    input:
        expand("sketches/metag/{n}.sig.zip", n=METAGENOME_NAMES),
        expand("sketches/genomes/{n}.sig.zip", n=GENOME_NAMES),


def genome_inp(wc):
    return GENOME_NAMES[wc.name]

rule sketch_genome:
    input:
        genome_inp
    output:
        "sketches/genomes/{name}.sig.zip"
    shell: """
        sourmash sketch dna {input} -o {output} \
           -p k=21,k=31,k=51,scaled=1000
    """


def metag_inp(wc):
    return METAGENOME_NAMES[wc.name]

rule sketch_metag:
    input:
        metag_inp
    output:
        "sketches/metag/{name}.sig.zip"
    shell: """
        sourmash sketch dna {input} -o {output} \
           -p abund,k=21,k=31,k=51,scaled=1000
    """
