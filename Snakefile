import glob, os

configfile: "config.yml"

wildcard_constraints:
    name='[^./]+',
    k = "\d+",

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

KSIZES = [21, 31, 51]
GATHER_KSIZE = 21
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
        expand("outputs/metag_compare.{k}.abund.matrix.png", k=KSIZES),
        expand("outputs/metag_compare.{k}.flat.matrix.png", k=KSIZES),
        expand("outputs/genome_compare.{k}.ani.matrix.png", k=KSIZES),
        expand("outputs/metag_gather/{n}.{k}.gather.txt",
               n=METAGENOME_NAMES, k=GATHER_KSIZE),
        expand("outputs/metag_gather/{n}.{k}.gather.csv",
               n=METAGENOME_NAMES, k=GATHER_KSIZE),

def genome_inp(wc):
    return GENOME_NAMES[wc.name]

rule sketch_genome:
    input:
        genome_inp
    output:
        "sketches/genomes/{name}.sig.zip"
    shell: """
        sourmash sketch dna {input:q} -o {output:q} \
           -p k=21,k=31,k=51,scaled=1000 \
           --name {name:q}
    """


def metag_inp(wc):
    return METAGENOME_NAMES[wc.name]

rule sketch_metag:
    input:
        metag_inp
    output:
        "sketches/metag/{name}.sig.zip"
    shell: """
        sourmash sketch dna {input:q} -o {output:q} \
           -p abund,k=21,k=31,k=51,scaled=1000 \
           --name {name:q}
    """

rule make_metagenome_compare_abund:
    input:
        expand("sketches/metag/{n}.sig.zip", n=METAGENOME_NAMES),
    output:
        cmp="outputs/metag_compare.{k}.abund",
        labels="outputs/metag_compare.{k}.abund.labels.txt"
    shell: """
        sourmash compare {input} -o {output.cmp} -k {wildcards.k}
    """

rule make_metagenome_compare_flat:
    input:
        expand("sketches/metag/{n}.sig.zip", n=METAGENOME_NAMES),
    output:
        cmp="outputs/metag_compare.{k}.flat",
        labels="outputs/metag_compare.{k}.flat.labels.txt"
    shell: """
        sourmash compare {input} -o {output.cmp} -k {wildcards.k} \
            --ignore-abund
    """

rule make_genome_compare_ani:
    input:
        expand("sketches/genomes/{n}.sig.zip", n=GENOME_NAMES),
    output:
        cmp="outputs/genome_compare.{k}.ani",
        labels="outputs/genome_compare.{k}.ani.labels.txt"
    shell: """
        sourmash compare {input} -o {output.cmp} -k {wildcards.k} \
            --containment --ani
    """

rule make_matrix_pdf:
    input:
        "outputs/{cmp}",
    output:
        "outputs/{cmp}.matrix.png"
    shell: """
        sourmash plot {input} --output-dir=outputs/
    """

rule unpack_database:
    input:
        config['databases']
    output:
        directory(config.get('database_sketches_dir', 'interim/database_sketches.d'))
    shell: """
        sourmash sig split {input} --outdir {output} -E .sig.gz \
            -k {GATHER_KSIZE}
    """
        
rule list_databases:
    input:
        config.get('database_sketches_dir', 'interim/database_sketches.d')
    output:
        "interim/list.database-sketches.txt"
    shell: """
        find {input} -name "*.sig.gz" > {output}
    """

rule metag_extract_sketch:
    input:
        "sketches/metag/{name}.sig.zip",
    output:
        "sketches/metag/{name}.{k}.sig.gz",
    shell: """
        sourmash sig cat {input} -o {output} -k {wildcards.k}
    """

rule metag_fastgather:
    input:
        query="sketches/metag/{name}.{k}.sig.gz",
        db="interim/list.database-sketches.txt"
    output:
        "outputs/metag_gather/{name}.{k}.fastgather.csv",
    shell: """
        sourmash scripts fastgather {input.query} {input.db} -o {output} \
           -k {wildcards.k}
    """

rule metag_gather:
    input:
        query = "sketches/metag/{name}.sig.zip",
        db = config['databases'],
        fastgather_out = "outputs/metag_gather/{name}.{k}.fastgather.csv",
    output:
        csv=touch("outputs/metag_gather/{name}.{k}.gather.csv"),
        out="outputs/metag_gather/{name}.{k}.gather.txt",
    shell: """
        sourmash gather {input.query} {input.db} -k {wildcards.k} \
            --picklist {input.fastgather_out}:match_md5:md5 \
            -o {output.csv} > {output.out}
    """
