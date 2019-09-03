import sys
import parse_checkm

import pysam

species = [
    "bacillus_subtilis",
    "enterococcus_faecalis",
    "escherichia_coli",
    "listeria_monocytogenes",
    "pseudomonas_aeruginosa",
    "saccharomyces_cerevisiae",
    "salmonella_enterica",
    "staphylococcus_aureus",
]

meta_fh = open(sys.argv[2])
fa_meta = {}
for line in meta_fh:
    fields = line.strip().split(",")
    fa_meta[fields[0]] = [fields[1], fields[2]]

manifest_fh = open(sys.argv[3])
uuid_header = manifest_fh.readline().strip().split("\t")
uuid_manifest = {}
for line in manifest_fh:
    fields = line.strip().split("\t")
    uuid_manifest[fields[0]] = dict(zip(uuid_header, fields))

try_checkm = False
if len(sys.argv) == 5:
    try_checkm = True

tab_fh = open("minidot_summary.tsv", 'w')

print("<table>")
fields = [
    "<tr>",
        "<th>Platform</th>",
        "<th>Zymo</th>",
        "<th>Extraction</th>",
        "<th>L</th>",
        "<th>e</th>",
        "<th>p</th>",

        "<th>Total bp</th>",
        "<th>Contigs</th>",
        "<th>N50</th>",
]

fields.extend([
    "<th>%s</th>" % "".join([x[0] for x in ref.split("_")]) for ref in species
])
fields.append("<th>Download</th>")

fields.append("</tr>")

print("\n".join(fields))

stats_fh = open(sys.argv[1])
for line in stats_fh:
    fields = line.strip().split('\t')
    if fields[0] == "Sample_ID":
        continue

    fa = fields[0]

    uuid = fa.split('.')[0]

    try:
        properties = uuid_manifest[uuid]
    except KeyError:
        continue

    genome_size = int(fields[1])
    n_contigs = int(fields[2])
    n50 = int(fields[5])

    fields = [
        "<tr>",
            "<td>%s</td>" % properties.get("platform", "?"),
            "<td>%s</td>" % properties.get("community", "?"),
            "<td>%s</td>" % properties.get("extraction", "?"),

            "<td>%s</td>" % properties.get("length", "?"),
            "<td>%s</td>" % properties.get("edge", "?"),
            "<td>%s</td>" % properties.get("pmer", "?"),

            "<td>%d</td>" % genome_size,
            "<td>%d</td>" % n_contigs,
            "<td>%d</td>" % n50,
    ]

    fields.extend([
        "<td><img height=150 width=150 src='minidotplots/%s.%s.png'></td>" % (fa.replace(".fa", ""), ref) for ref in species
    ])

    fields.append("<td><a href='http://nanopore.s3.climb.ac.uk/mockcommunity/v3/%s'>Download (%s, <code>%s</code>)</a></br>%s</td>" % (fa, fa_meta[fa][0], fa_meta[fa][1], fa.replace(uuid, "")))
    fields.append("</tr>")

    if try_checkm:
        fields.append("<tr><td colspan=9></td>")

        checkm_d = {}
        checkm_d = parse_checkm.parse_checkm("checkm-%s" % fa.replace(".fa", ".txt"))
        fields.extend([
            "<td style='text-align:center'>%s</td>" % (checkm_d.get(ref, {}).get("Completeness", "?")) for ref in species
        ])
        fields.append("</tr>")

    for ref in species:
        try:
            # To be honest we could have just done this from the FAI but whatever
            extracted_fa = pysam.FastaFile("extracted_contigs/%s/%s.fasta" % (fa.replace(".fa", ""), ref))

            len_list = extracted_fa.lengths
            genome_total_size = sum(len_list)
            genome_largest_contig = max(len_list)
            genome_n_contigs = len(len_list)

            extracted_fa.close()
        except:
            genome_total_size = 0
            genome_largest_contig = 0
            genome_n_contigs = 0

        tab_fh.write("\t".join([ str(x) for x in [
            uuid,
            fa.replace(uuid, ""),
            ref,
            genome_total_size,
            genome_largest_contig,
            genome_n_contigs,
            checkm_d.get(ref, {}).get("Completeness", "?"),
        ]]) + '\n')

    print("\n".join(fields))

print("</table>")

tab_fh.close()
