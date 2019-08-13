import sys

gfa = open(sys.argv[1])
LIMIT = int(sys.argv[2])
#filt = open(sys.argv[3])

contigs_not_interest = {}
#for line in filt:
#    fields = line.strip().split()
#    name = fields[0].replace("contig", "edge")
#    contigs_not_interest[name] = 1

contigs_in_size = []

for line in gfa:
    fields = line.strip().split()

    if line[0] == 'S':
        name = fields[1]
        curr_len = None
        if fields[2] == "*":
            for meta in fields[3].split():
                if meta.startswith("LN"):
                    curr_len = int(meta.split(":")[2])
        else:
            curr_len = len(fields[2])
            name = fields[1]

        if curr_len >= LIMIT and name not in contigs_not_interest:
            contigs_in_size.append(name)
            print(line.strip())

    if line[0] == 'L':
        left = fields[1]
        right = fields[3]
        if left in contigs_in_size and right in contigs_in_size:
            print(line.strip())

