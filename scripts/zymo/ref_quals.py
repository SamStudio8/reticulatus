"""Given an indexed FASTQ and a BAM of aligned reads, output their name, ref, size and qual

Usage: <fai> <bam>"""
import sys
import pysam

#from nanomath import ave_qual

uuids = set([])

def convert_qualstr_to_ints(s):
    return [ord(c)-33 for c in s]

fh = open( sys.argv[1] )
for line in fh:
    fields = line.strip().split("\t")
    uuids.add(fields[0])

sys.stderr.write("completed %d uuids\n" % len(uuids))
data = dict.fromkeys(uuids) # presize dict

fh.seek(0)
for line in fh:
    fields = line.strip().split("\t")
    data[fields[0]] = {
        "len": int(fields[1]),
        "seek_seq": int(fields[2]),
        "seek_qual": int(fields[5]),
    }

sys.stderr.write("completed %d uuids\n" % len(data))
fq = open( sys.argv[1].replace(".fai", "") )


bam_fh = pysam.AlignmentFile( sys.argv[2] )

print("uuid ref len qual rlen")
for ref in bam_fh.references:
    ref_code = ref.split("_")[0]
    seen_on_ref = set([])

    for read in bam_fh.fetch(contig=ref):

        if read.mapping_quality == 0:
            continue

        if read.query_name in seen_on_ref:
            continue

        d = read.query_name
        #:fq.seek(data[d]["seek_qual"])
        print (d, ref_code, data[d]["len"], -1, read.reference_length) #ave_qual( convert_qualstr_to_ints(fq.read(data[d]["len"]))) )
        seen_on_ref.add(read.query_name)


fh.close()
fq.close()
