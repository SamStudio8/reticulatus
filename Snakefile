import sys

# Import config, set working dir, load manifest and define default rules:
include: "Snakefile-base"

GENOMES = [
    "bacillus_subtilis",
    "enterococcus_faecalis",
    "escherichia_coli",
    "listeria_monocytogenes",
    "pseudomonas_aeruginosa",
    "saccharomyces_cerevisiae",
    "salmonella_enterica",
    "staphylococcus_aureus",
]

rule all:
    input:
        web_html="index.html",
        web_pngs="png.tar",
        kraken="kraken_summary.bond.tsv",
        polished=enumerate_assemblies(unroll=False, suffix=".fa")

rule mkdir_checkm:
    output: directory('checkm')
    shell: 'mkdir checkm; cd checkm; wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz; tar xvf checkm_data_2015_01_16.tar.gz'

rule setup_checkm:
    input: directory('checkm')
    conda: "environments/checkm.yaml"
    output: "checkm_setup.ok"
    shell: "python --version; echo 'checkm/' | checkm data setRoot 'checkm/'; touch checkm_setup.ok;"

def checkm_pick_taxon(w):
    lookup = {
        "bacillus_subtilis":        {"rank": "species", "genome": "Bacillus subtilis"},
        "enterococcus_faecalis":    {"rank": "species", "genome": "Enterococcus faecalis"},
        "escherichia_coli":         {"rank": "species", "genome": "Escherichia coli"},
        "listeria_monocytogenes":   {"rank": "species", "genome": "Listeria monocytogenes"},
        "pseudomonas_aeruginosa":   {"rank": "species", "genome": "Pseudomonas aeruginosa"},
        #"saccharomyces_cerevisiae": {"rank": "species", "genome": ""},
        "salmonella_enterica":      {"rank": "species", "genome": "Salmonella enterica"},
        "staphylococcus_aureus":    {"rank": "species", "genome": "Staphylococcus aureus"},
    }
    try:
        return lookup[w.genome]
    except KeyError:
        return {"rank": "life", "genome": "Prokaryote"}


# This rule is a little annoying, we *should* be explicitly naming the input
# FASTA needed for the checkM bin as extracted by our kraken script. Unfortunately,
# as the extract_kraken_contigs rule output is defined as a dictionary, we cannot
# explicitly name an input inside that directory without a ChildIOException.
# That rule has a python script that attempts to guarantee those files exist, so
# as long as the directory exists, we can be reasonably confident the FASTA do too.
#   As gross as it sounds, the requisite `cp` will fail if it *really* doesn't exist.
# So this isn't a particularly unsafe workaround, just a gross one.
rule checkm:
    input:
        ok="checkm_setup.ok",
        d=directory("extracted_contigs/{assembly}/")
        #fa="extracted_contigs/{assembly}/{genome}.fasta"
    output:
        working=directory("checkm-working/{assembly}/{genome}"),
        bin=directory("checkm-bins/{assembly}/{genome}/"),
        log="checkm-results/{assembly}/{genome}.tsv",
    params: p=checkm_pick_taxon
    conda: "environments/checkm.yaml"
    threads: 8
    shell: "cp extracted_contigs/{wildcards.assembly}/{wildcards.genome}.fasta {output.bin}; checkm taxonomy_wf -t {threads} -x fasta {params.p[rank]} {params.p[genome]:q} {output.bin} {output.working} > {output.log}"

rule merge_checkm:
    input:
        logs=expand("checkm-results/{{assembly}}/{genome}.tsv", genome=GENOMES)
    output:
        "checkm-{assembly}.txt"
    shell:
        "cat {input.logs} > {output}"

rule tabulate_minidot:
    input:
        stat="assembly_stats.txt",
        meta="assembly_md5size.txt",
        manifest="../manifest.cfg",
        dots=expand("minidotplots/{assembly}.{genome}.png", assembly=enumerate_assemblies(), genome=GENOMES),
        checkm=expand("checkm-{assembly}.txt", assembly=enumerate_assemblies())
    output:
        html="index.html",
        pngs="png.tar"
    shell:
        "python ../scripts/summarise_assemblies.py {input.stat} {input.meta} {input.manifest} 1 > {output.html}; tar -cvf {output.pngs} minidotplots/*.png"

rule install_minidot:
    output: touch("minidot.ok")
    shell: "cd git; git clone https://github.com/samstudio8/minidot.git"

rule minidot:
    input:
        d=directory("extracted_contigs/{assembly}/"),
        ref="mason_refs/{genome}_pb.fasta",
        ready="minidot.ok"
    output:
        "minidotplots/{assembly}.{genome}.png"
    shell:
        "python git/minidot/minidot.py --ignore-missing --mapper minimap2 --no-self --alen 2500 --strip git/minidot/minidot.R {output} {input.d}{wildcards.genome}.fasta {input.ref}"

rule extract_kraken_contigs:
    input:
        assembly="{assembly}.fa",
        k2report="{assembly}.fa.k2",
    params:
        genomes=GENOMES,
    output:
        directory("extracted_contigs/{assembly}/")
    shell:
        "python ../scripts/extract_contigs_with_kraken.py {input.k2report} {input.assembly} {output}; python ../scripts/ensure_genomes.py {output} {params.genomes}"

rule summarise_kraken:
    input:
        enumerate_assemblies(suffix=".fa.k2"), "static/Zymo-Isolates-SPAdes-Illumina.fasta.k2", "static/expected_genomes.k2", "static/pacbio.k2"
    output:
        "kraken_summary.tsv"
    shell:
        "python ../scripts/extracken2.py {input} > {output}"

