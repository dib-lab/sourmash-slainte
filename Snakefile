import glob, os, csv
from collections import defaultdict

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
    return basename

KSIZES = [21, 31, 51]
GATHER_KSIZE = 21

GENOME_NAMES = {}
for g in config['genomes']:
    files = glob.glob(g)
    for filename in files:
        name = strip_suffix(filename)
        assert name not in GENOME_NAMES, f"duplicate prefix for {filename}"
        GENOME_NAMES[name] = filename

print(f"Found {len(GENOME_NAMES)} genome files.")
if len(GENOME_NAMES) > 0:
    ENABLE_GENOMES = True
else:
    print('** NOTE: no genome files found. Disabling genome output!')
    ENABLE_GENOMES = False

METAG_PATH=config['metagenome_dir']
METAGENOME_NAMES=defaultdict(set)
with open(config['sample_info'], 'r', newline='') as sample_fp:
    r = csv.DictReader(sample_fp)

    for row in r:
        fileglob = METAG_PATH.rstrip('/') + '/' + row['prefix'] + '*'
        name = row['name']

        files = glob.glob(fileglob)
        print('F', name, fileglob, files)
        assert files, fileglob

        METAGENOME_NAMES[name].update(files)

print(f"Found {len(METAGENOME_NAMES)} metagenome names.")


genome_outputs = []
if ENABLE_GENOMES:
    genome_outputs.extend(
        expand("outputs/genome_compare.{k}.ani.matrix.png", k=KSIZES)
        )
    genome_outputs.extend(
        expand("outputs/metag.x.genomes.{k}.png", k=KSIZES),
        )
    genome_outputs.extend(
        expand("outputs/prefetch/all_metag.x.genomes.{k}.summary.png",
               k=GATHER_KSIZE),
        )

rule all:
    input:
        expand("sketches/metag/{n}.sig.zip", n=METAGENOME_NAMES),
        expand("sketches/genomes/{n}.sig.zip", n=GENOME_NAMES),
        expand("outputs/metag_compare.{k}.abund.matrix.png", k=KSIZES),
        expand("outputs/metag_compare.{k}.flat.matrix.png", k=KSIZES),
        expand("outputs/metag_gather/{n}.{k}.gather.txt",
               n=METAGENOME_NAMES, k=GATHER_KSIZE),
        expand("outputs/metag_gather/{n}.{k}.gather.csv",
               n=METAGENOME_NAMES, k=GATHER_KSIZE),
        expand("outputs/metag_gather/{n}.{k}.kreport.out",
               n=METAGENOME_NAMES, k=GATHER_KSIZE),
        expand("outputs/metag_gather/metag.{k}.kreport.csv", k=GATHER_KSIZE),
        genome_outputs,


rule sketch:
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
        sourmash sketch dna {input:q} -o {output:q} \
           -p k=21,k=31,k=51,scaled=1000 \
           --name {wildcards.name:q}
    """


def metag_inp(wc):
    return list(METAGENOME_NAMES[wc.name])

rule sketch_metag:
    input:
        metag_inp
    output:
        "sketches/metag/{name}.sig.zip"
    shell: """
        sourmash sketch dna {input:q} -o {output:q} \
           -p abund,k=21,k=31,k=51,scaled=1000 \
           --name {wildcards.name:q}
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

rule extract_sketch:
    input:
        "sketches/{dir}/{name}.sig.zip",
    output:
        "sketches/{dir}/{name}.{k}.sig.gz",
    shell: """
        sourmash sig cat {input} -o {output} -k {wildcards.k}
    """

rule metag_fastgather:
    input:
        query="sketches/metag/{name}.{k}.sig.gz",
        db="interim/list.database-sketches.txt"
    output:
        "outputs/metag_gather/{name}.{k}.fastgather.csv",
    threads: 64
    shell: """
        sourmash scripts fastgather {input.query} {input.db} -o {output} \
           -k {wildcards.k} -c {threads}
    """

rule metag_gather:
    input:
        query = "sketches/metag/{name}.sig.zip",
        db = config['databases'],
        fastgather_out = "outputs/metag_gather/{name}.{k}.fastgather.csv",
    output:
        csv=touch("outputs/metag_gather/{name}.{k}.gather.csv"),
        out="outputs/metag_gather/{name}.{k}.gather.txt",
    resources:
        gather=1
    shell: """
        sourmash gather {input.query} {input.db} -k {wildcards.k} \
            --picklist {input.fastgather_out}:match_md5:md5 \
            -o {output.csv} > {output.out}
    """

rule prepare_taxdb:
    input:
        config['taxonomies']
    output:
        "interim/tax.sqldb"
    shell: """
        sourmash tax prepare -t {input} -o {output} -F sql
    """

rule metag_tax:
    input:
        gather_csv = "outputs/metag_gather/{name}.{k}.gather.csv",
        taxdb = "interim/tax.sqldb",
    output:
        "outputs/metag_gather/{name}.{k}.kreport.out",
    shell: """
        sourmash tax metagenome -F kreport -g {input.gather_csv} \
            -t {input.taxdb} > {output}
    """

rule metag_tax_summary:
    input:
        expand("outputs/metag_gather/{n}.{{k}}.kreport.out",
               n=METAGENOME_NAMES),
    output:
        "outputs/metag_gather/metag.{k}.kreport.csv"
    shell: """
        scripts/combine-kreports.py {input} -o {output}
    """

rule list_genomes:
    input:
        expand("sketches/genomes/{n}.{{k}}.sig.gz", n=GENOME_NAMES),
    output:
        "interim/list.genome-sketches.{k}.txt"
    run:
        with open(output[0], 'wt') as outfp:
            outfp.write("\n".join(input))

rule list_metag:
    input:
        expand("sketches/metag/{n}.{{k}}.sig.gz", n=METAGENOME_NAMES),
    output:
        "interim/list.metag-sketches.{k}.txt"
    shell: """
        ls -1 {input} > {output}
    """

rule metag_x_genomes_csv:
    input:
        metag="interim/list.metag-sketches.{k}.txt",
        genomes="interim/list.genome-sketches.{k}.txt",
    output:
        "outputs/metag.x.genomes.{k}.manysearch.csv"
    threads: 8
    shell: """
        sourmash scripts manysearch -k {wildcards.k} \
            {input.genomes} {input.metag} \
            -c {threads} -t 0 \
            -o {output}
    """

rule summarize_manysearch:
    input:
        "outputs/metag.x.genomes.{k}.manysearch.csv"
    output:
        "outputs/metag.x.genomes.{k}.summary.csv"
    shell: """
        scripts/summarize-manysearch.py {input} -o {output}
    """

rule plot_manysearch:
    input:
        "outputs/metag.x.genomes.{k}.summary.csv"
    output:
        "outputs/metag.x.genomes.{k}.png"
    shell: """
        scripts/plot-genome-vs-metag.py {input} -o {output}
    """

rule metag_x_genomes_prefetch:
    input:
        genomes = expand("sketches/genomes/{n}.sig.zip", n=GENOME_NAMES),
        metag="sketches/metag/{metag}.sig.zip",
        bin = "scripts/calc-weighted-overlap.py",
    output:
        "outputs/prefetch/{metag}.x.genomes.{k}.csv",
    threads: 1
    shell: """
        {input.bin} -k {wildcards.k} --genomes {input.genomes} \
            --metagenomes {input.metag} -o {output}
    """

rule metag_x_genomes_prefetch_summary:
    input:
        csv=expand("outputs/prefetch/{metag}.x.genomes.{k}.csv",
               metag=METAGENOME_NAMES, k=GATHER_KSIZE),
        bin="scripts/summarize-weighted-overlap.py",
    output:
        "outputs/prefetch/all_metag.x.genomes.{k}.summary.csv",
    shell: """
        {input.bin} {input.csv} -o {output}
    """

rule plot_prefetch_summary:
    input:
        "outputs/prefetch/all_metag.x.genomes.{k}.summary.csv"
    output:
        "outputs/prefetch/all_metag.x.genomes.{k}.summary.png"
    shell: """
        scripts/plot-genome-vs-metag.py {input} -o {output}
    """
