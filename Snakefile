# import utility functions:
from slainte_functions import *

configfile: "config.yml"

wildcard_constraints:
    name='[^./]+',              # for 'name', don't recurse into subdirectories
    k = "\\d+",                 # for k-mer sizes, allow numbers only

# k-mer sizes to sketch:
KSIZES = [21, 31, 51]

# param strings for sketching
k_params = ",".join(expand("k={k}", k=KSIZES))
genome_sketch_params = "scaled=1000," + k_params
metag_sketch_params = "scaled=1000,abund," + k_params

# k-mer sizes for running gather:
GATHER_KSIZE = 21

#
# collect all genome files.
#

GENOME_NAMES = collect_genomes(config['genomes'])

if len(GENOME_NAMES) > 0 and config['display_genomes']:
    ENABLE_GENOMES = True
else:
    print('** NOTE: no genome files found OR display_genomes is False. Disabling genome output!')
    ENABLE_GENOMES = False

#
# collect all metagenome files, based on contents of sample_info CSV file.
#

metag_path = config['metagenome_dir']
samples_csv = config['sample_info']
METAGENOME_NAMES, METAGENOME_FILES = load_metagenome_files(metag_path,
                                                           samples_csv,
                                                           debug=True)

print(f"Found {len(METAGENOME_NAMES)} samples total!")

#
# configure - run gather?
#

if config['run_gather']:
    RUN_GATHER = True
else:
    print('** NOTE: run_gather is False. Disabling gather output!')
    RUN_GATHER = False

#####
##### all configuration decisions made! now, make it so.
#####

extra_outputs = []
if ENABLE_GENOMES:
    extra_outputs.extend(
        expand("outputs/genome_compare.{k}.ani.matrix.png", k=KSIZES)
        )
    extra_outputs.extend(
        expand("outputs/metag.x.genomes.{k}.manysearch.png", k=KSIZES),
        )
    extra_outputs.extend(
        expand("outputs/prefetch/all_metag.x.genomes.{k}.summary.png",
               k=GATHER_KSIZE),
        )

if RUN_GATHER:
    extra_outputs.extend(
        expand("outputs/metag_gather/{n}.{k}.gather.csv",
               n=METAGENOME_NAMES, k=GATHER_KSIZE),
        )
    extra_outputs.extend(
        expand("outputs/metag_gather/metag.{k}.kreport.csv", k=GATHER_KSIZE),
        )

rule all:
    input:
        expand("sketches/metag/{n}.sig.zip", n=METAGENOME_NAMES),
        expand("sketches/metag_individual_sketches/{f}.sig.zip",
               f=METAGENOME_FILES),
        "outputs/check/metag.x.individual.21.manysearch.png",
        expand("sketches/genomes/{n}.sig.zip", n=GENOME_NAMES),
        expand("outputs/metag_compare.{k}.abund.matrix.png", k=KSIZES),
        expand("outputs/metag_compare.{k}.flat.matrix.png", k=KSIZES),
        extra_outputs,


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
           -p {genome_sketch_params} --name {wildcards.name:q}
    """

def metag_individual_inp(wc):
    return METAGENOME_FILES[wc.name]

rule sketch_metag_individual_data_file:
    input:
        metag_individual_inp
    output:
        "sketches/metag_individual_sketches/{name}.sig.zip",
    shell: """
        echo name,genome_filename,protein_filename > {output:q}.manysketch.csv
        echo {wildcards.name},{input:q}, >> {output:q}.manysketch.csv
        sourmash scripts manysketch {output:q}.manysketch.csv -o {output:q} \
           -p {metag_sketch_params} -c 1

    """

def metag_get_individual_sketches(wc):
    name = wc.name              # this will be sample name for combined file

    # now get the individual sketch files for this sample
    datafiles = METAGENOME_NAMES[name]
    sketch_names = []
    for datafile in datafiles:
        for k, v in METAGENOME_FILES.items():
            if datafile == v:
                sketch_names.append(k)

    # convert to filenames
    return expand("sketches/metag_individual_sketches/{f}.sig.zip",
                  f=sketch_names)

rule merge_individual_sketches:
    input:
        metag_get_individual_sketches
    output:
        "sketches/metag/{name}.{k}.sig.gz"
    shell: """
        sourmash sig merge -k {wildcards.k} {input:q} -o {output:q} \
           --name {wildcards.name:q}
    """

rule combine_merged_sketches:
    input:
        expand("sketches/metag/{{name}}.{k}.sig.gz", k=KSIZES)
    output:
        "sketches/metag/{name}.sig.zip"
    shell: """
        sourmash sig cat {input} -o {output}
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

rule extract_genome_sketch:
    input:
        "sketches/genomes/{name}.sig.zip",
    output:
        "sketches/genomes/{name}.{k}.sig.gz",
    shell: """
        sourmash sig cat {input} -o {output} -k {wildcards.k}
    """

rule extract_individual_sketch:
    input:
        "sketches/metag_individual_sketches/{name}.sig.zip",
    output:
        "sketches/metag_individual_sketches/{name}.{k}.sig.gz",
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

rule list_individual_metag_sketches:
    input:
        expand("sketches/metag_individual_sketches/{n}.{{k}}.sig.gz", n=METAGENOME_FILES),
    output:
        "interim/list.metag-individual-sketches.{k}.txt"
    shell: """
        ls -1 {input} > {output}
    """

rule metag_x_individual_sketches_csv:
    input:
        metag="interim/list.metag-sketches.{k}.txt",
        individual="interim/list.metag-individual-sketches.{k}.txt",
    output:
        "outputs/check/metag.x.individual.{k}.manysearch.raw.csv"
    threads: 8
    shell: """
        sourmash scripts manysearch -k {wildcards.k} \
            {input.individual} {input.metag} \
            -c {threads} -t 0 \
            -o {output}
    """

rule metag_x_genomes_csv:
    input:
        metag="interim/list.metag-sketches.{k}.txt",
        genomes="interim/list.genome-sketches.{k}.txt",
    output:
        "outputs/metag.x.genomes.{k}.manysearch.raw.csv"
    threads: 8
    shell: """
        sourmash scripts manysearch -k {wildcards.k} \
            {input.genomes} {input.metag} \
            -c {threads} -t 0 \
            -o {output}
    """

rule summarize_manysearch:
    input:
        "{path}.manysearch.raw.csv"
    output:
        "{path}.manysearch.summary.csv"
    shell: """
        scripts/summarize-manysearch.py {input} -o {output}
    """

rule plot_manysearch:
    input:
        "{path}.manysearch.summary.csv"
    output:
        "{path}.manysearch.png"
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
