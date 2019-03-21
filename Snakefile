import pandas as pd
from snakemake.utils import validate
import sys

shell.executable('bash')

# Import config
configfile: "config.yaml"

# Set working dir to "working/", you can change this to anything.
# The rationale is you can git clone the pipeline and run it without doing any juggling.
# Note that most things need to be prefixed with ../ to refer to the current directory,
# but there are exceptions like when reading cfg files and conda environments.
workdir: "working/"

# Load manifest(s)
samples = pd.read_table("../manifest.cfg").set_index("uuid", drop=False)     # defines the assembly and polishing strategies
illumina_lookup = pd.read_table("../reads.cfg").set_index("ont", drop=False) # provides a lookup of short reads for polishing

# Utility functions
def unroll_assemblies(w, name, unroll=True):
    if not unroll:
        return [name]

    polishes = name.split(".")
    if len(polishes) == 0 or len(polishes[0]) == 0:
        return []

    end_polish = polishes[-1]
    end_polish_f = end_polish.split("-")
    end_iter = int(end_polish_f[-1])
    if (end_iter - 1) == 0:
        return [name] + unroll_assemblies(w, ".".join(polishes[:-1]))
    else:
        end_polish_f[-1] = str(end_iter-1)
        return [name] + unroll_assemblies(w, ".".join(polishes[:-1] + ["-".join(end_polish_f)]))

def enumerate_assemblies(w=None, suffix="", unroll=True):
    base_assemblies = []
    polished_assemblies = []
    for uuid in samples["uuid"]:
        for assembler in samples.loc[uuid]["assemblers"].split(","):
            base_assemblies.append( "%s.%s.ctg.cns%s" % (uuid, assembler, suffix) )
            if samples.loc[uuid]["poa"] != '-':
                unrolled = unroll_assemblies(w, samples.loc[uuid]["poa"], unroll=unroll)
                context = ['%s.%s.ctg.cns.%s%s' % (uuid, assembler, step, suffix) for step in unrolled]
                polished_assemblies.extend(context)
    return base_assemblies + polished_assemblies

for a in enumerate_assemblies(suffix=".fa"):
    sys.stderr.write("*\t%s\n" % a)

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

rule finish:
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

#rule get_ready_to_Rumble:
#    input:
#    output:
#        ""
#    shell:
#        "Rscript  --vanilla"

def input_polish(w, name):
    if int(w.iteration) > 1:
        contigs = "%s.%s.ctg.cns.%s%s-%s-%d.fa" % (w.uuid, w.assembler, w.polishedprefix, name, w.readtype, int(w.iteration)-1)
    elif w.polishedprefix=="":
        contigs = "%s.%s.ctg.cns.fa" % (w.uuid, w.assembler)
    else:
        contigs = "%s.%s.ctg.cns.%sfa" % (w.uuid, w.assembler, w.polishedprefix)
    return contigs

def input_polish_racon(w):
    return input_polish(w, "racon")
def input_polish_medaka(w):
    return input_polish(w, "medaka")
def input_polish_pilon(w):
    return input_polish(w, "pilon")

def minimap2_mode(w):
    if w.readtype == "ont":
        mode = "map-ont"
    elif w.readtype == "ill":
        mode = "sr"
    return mode

def polish_reads_input(w):
    if w.readtype == "ont":
        reads = os.path.join(config["long_fq_root"], samples.loc[w.uuid]['reads'])
    elif w.readtype == "ill":
        reads = expand(os.path.join(config["short_fq_root"], "{fq}"), fq=[illumina_lookup.loc[ samples.loc[w.uuid]['reads'] ]["i1"], illumina_lookup.loc[ samples.loc[w.uuid]['reads'] ]["i2"]])
    return reads

rule polish_racon:
    input: contigs=input_polish_racon, reads=polish_reads_input
    params:
        mode=minimap2_mode,
    output:
        "{uuid}.{assembler}.ctg.cns.{polishedprefix,.*}racon-{readtype,\w+}-{iteration,\d+}.fa"
    threads: 12
    shell:
        "minimap2 -t {threads} -x {params.mode} {input.contigs} {input.reads} > {output}.paf; racon -t {threads} {input.reads} {output}.paf {input.contigs} > {output}"

rule polish_medaka:
    input: contigs=input_polish_medaka, reads=polish_reads_input
    params:
        model=lambda w: samples.loc[w.uuid]['medakamodel'],
    output:
        "{uuid}.{assembler}.ctg.cns.{polishedprefix,.*}medaka-{readtype,\w+}-{iteration,\d+}.fa"
    threads: 12
    shell:
        "rm -rf medaka-{wildcards.uuid}/*; medaka_consensus -i {input.reads} -d {input.contigs} -o medaka-{wildcards.uuid} -t {threads} -m {params.model}; mv medaka-{wildcards.uuid}/consensus.fasta {output}"

rule download_pilon:
    output: "pilon-1.23.jar"
    shell: "wget https://github.com/broadinstitute/pilon/releases/download/v1.23/pilon-1.23.jar"

rule polish_pilon:
    input: contigs=input_polish_pilon, reads=polish_reads_input, pilon="pilon-1.23.jar"
    params:
        mode=minimap2_mode,
    output:
        polish="{uuid}.{assembler}.ctg.cns.{polishedprefix,.*}pilon-{readtype,\w+}-{iteration,\d+}.fa"
    threads: 12
    shell:
        "echo {input.reads}; minimap2 -t {threads} -ax {params.mode} {input.contigs} {input.reads} > {output}.sam; samtools sort {output}.sam -T {wildcards.uuid} -m 2G -@ {threads} -o {output}.bam; samtools index {output}.bam; java -Xmx16G -jar {input.pilon} --genome {input.contigs} --bam {output}.bam --outdir pilon-{wildcards.uuid}/; mv pilon-{wildcards.uuid}/pilon.fasta {output}"

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
        "python scripts/summarise_assemblies.py {input.stat} {input.meta} {input.manifest} 1 > {output.html}; tar -cvf {output.pngs} minidotplots/*.png"

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
        "python scripts/extract_contigs_with_kraken.py {input.k2report} {input.assembly} {output}; python scripts/ensure_genomes.py {output} {params.genomes}"

rule summarise_assembly_stats:
    input:
        enumerate_assemblies(suffix=".fa")
    output:
        "assembly_stats.txt"
    shell:
        "for fa in {input}; do perl scripts/assembly-stats.pl $fa; done > {output}"

rule summarise_assembly_meta:
    input:
        enumerate_assemblies(suffix=".fa")
    output:
        "assembly_md5size.txt"
    shell:
        "for fa in {input}; do base=`basename $fa`; md5=`md5sum $fa | cut -f1 -d' '`; size=`du -h $fa | cut -f1`; echo \"$base,$size,$md5\"; done > {output}"

rule bond_summarise_kraken:
    input:
        "kraken_summary.tsv"
    output:
        "kraken_summary.bond.tsv"
    shell:
        "bond {input} samples.cfg > {output}"

rule summarise_kraken:
    input:
        enumerate_assemblies(suffix=".fa.k2"), "static/Zymo-Isolates-SPAdes-Illumina.fasta.k2", "static/expected_genomes.k2", "static/pacbio.k2"
    output:
        "kraken_summary.tsv"
    shell:
        "python scripts/extracken2.py {input} > {output}"

rule download_kraken_database:
    output:
        ok=touch("k2db.ok"),
        h=os.path.join(config["database_root"], "hash.k2d"),
        o=os.path.join(config["database_root"], "opts.k2d"),
        t=os.path.join(config["database_root"], "taxo.k2d"),
    shell: "cd %s; wget https://raw.githubusercontent.com/LomanLab/mockcommunity/master/analysis/snakemake-assembly/databases/go.sh; bash go.sh" % config["database_root"]

rule kraken:
    input:
        fa="{uuid}.{prefix}.fa", ok="k2db.ok"
    output:
        "{uuid}.{prefix}.fa.k2"
    threads: 8
    shell:
        "kraken2 --db %s --use-names -t {threads} --output {output} {input.fa}" % config["database_root"]

def pick_assembler_version(assembler):
    lookup = {
        "wtdbg2": "bin/wtdbg2",
        "wtdbg2_2-4": "bin/wtdbg2_2-4",
    }
    return lookup[assembler]

def pick_wtdbg2_cns_version(assembler):
    lookup = {
        "wtdbg2": "bin/wtpoa-cns",
        "wtdbg2_2-4": "bin/wtpoa-cns_2-4",
    }
    return lookup[assembler]

rule wtdbg2_consensus:
    input:
        "{uuid}.{assembler}.ctg.lay.gz"
    output:
        "{uuid}.{assembler,wtdbg2[A-z0-9_-]*}.ctg.cns.fa"
    params:
        cnsbin=lambda w: pick_wtdbg2_cns_version(w.assembler)
    threads: 12
    shell:
        "{params.cnsbin} -f -i {input} -o {output} -t {threads}"

# new and improved version of this rule doesn't fucking nuke your input reads
rule subset_reads:
    input:
        os.path.join(config["long_fq_root"], "{fq}.fq.gz")
    output:
        os.path.join(config["long_fq_root"], "{fq}.subset.{ratio}.fq.gz")
    params:
        ratio=lambda w: float(w.ratio)/100
    shell:
        "seqtk sample {input} {params.ratio} | gzip -1 > {output}"

rule install_wtdbg2:
    output:
        ok=touch("w2.ok"),
        d=directory("git/wtdbg2"),
        w2_bin_904f2b3="bin/wtdbg2",
        cns_bin_904f2b3="bin/wtpoa-cns",
        w2_bin_6a0691e="bin/wtdbg2_2-4",
        cns_bin_6a0691e="bin/wtpoa-cns_2-4"
    shell: "cd git; git clone https://github.com/ruanjue/wtdbg2.git; cd wtdbg2; git reset --hard 904f2b3ebdaa1e6f268cc58937767891a00d5bcb; make; cp wtdbg2 ../../bin; cp wtpoa-cns ../../bin; git reset --hard 6a0691e308b3644b6f718a03679f697d058e2be6; make; cp wtdbg2 ../../bin/wtdbg2_2-4; cp wtpoa-cns ../../bin/wtpoa-cns_2-4;"

rule wtdbg2_assembly:
    input:
        reads=lambda w: os.path.join(config["long_fq_root"], samples.loc[w.uuid]['reads']),
        ready="w2.ok"
    params:
        abin=lambda w: pick_assembler_version(w.assembler),
        prefix=lambda w: samples.loc[w.uuid]['uuid'],
        pmer=lambda w: samples.loc[w.uuid]['pmer'],
        sampler=lambda w: samples.loc[w.uuid]['sampler'],
        edge=lambda w: samples.loc[w.uuid]['edge'],
        length=lambda w: samples.loc[w.uuid]['length'],
        max_k=lambda w: samples.loc[w.uuid]['max_k'],
        max_node=lambda w: samples.loc[w.uuid]['max_node'],
    output:
        "{uuid}.{assembler,wtdbg2[A-z0-9_-]*}.ctg.lay.gz"
    threads: 24
    shell:
        "{params.abin} -f -i {input.reads} -o {params.prefix} -S {params.sampler} -e {params.edge} -p {params.pmer} -L {params.length} -K {params.max_k} --node-max {params.max_node} -t {threads}"

rule install_flye:
    output:
        ok=touch("flye.ok"),
        d=directory("git/Flye"),
        flye_bin="bin/flye",
    shell: "cd git; git clone https://github.com/fenderglass/Flye.git; cd Flye; git reset --hard b3cfc9d96798307dfe91bfc1543de0dcba944cbf; make; cp bin/flye ../../bin;"

# TODO We need to use git/Flye/bin/flye as the script sets up some dir dependent stuff
rule flye_assembly:
    conda: "environments/flye.yaml"
    input:
        reads=lambda w: os.path.join(config["long_fq_root"], samples.loc[w.uuid]['reads']),
        ready="flye.ok"
    params:
        genomesize="62m",
        d = "{uuid}.{assembler,flye[A-z0-9_-]*}/"
    output:
        fa = "{uuid}.{assembler,flye[A-z0-9_-]*}/assembly.fasta",
        gfa = "{uuid}.{assembler,flye[A-z0-9_-]*}/assembly.gfa"
    threads: 48
    shell:
        "git/Flye/bin/flye --nano-raw {input.reads} --meta --plasmids -g {params.genomesize} -o {params.d} -t {threads}"

rule link_flye_assembly:
    input: "{uuid}.{assembler}/assembly.fasta"
    output: "{uuid}.{assembler,flye[A-z0-9._-]*}.ctg.cns.fa"
    shell: "ln -s {input} {output}"
