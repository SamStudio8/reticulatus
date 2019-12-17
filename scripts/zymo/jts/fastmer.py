#! /usr/bin/env python

import os
import pysam
import argparse
import vcf
import math
import sys
import subprocess
import numpy as np
from collections import namedtuple
from collections import defaultdict

# class to summarize a pairwise alignment between two sequences
# like a genome assembly to a reference genome
class AlignmentStats():
    
    def __init__(self):
        self.matches = 0
        self.mismatches = 0
        self.insertions = 0
        self.deletions = 0
        self.hp_count = defaultdict(int)
        self.hp_correct = defaultdict(int)

    def add_match(self):
        self.matches += 1
    
    def add_mismatch(self):
        self.mismatches += 1
    
    def add_insertion(self):
        self.insertions += 1

    def add_deletion(self):
        self.deletions += 1

    def get_identity(self):
        alignment_length = self.matches + self.mismatches + self.insertions + self.deletions
        return float(self.matches) / alignment_length

class AssemblyAccuracy():
    
    def __init__(self):
        self.window_size = 100000
        self.stats_by_window = defaultdict(AlignmentStats)
        self.global_stats = AlignmentStats()

    def add_match(self, position):
        self.global_stats.add_match()
        self.stats_by_window[self.get_bin(position)].add_match()
    
    def add_mismatch(self, position):
        self.global_stats.add_mismatch()
        self.stats_by_window[self.get_bin(position)].add_mismatch()
    
    def add_insertion(self, position):
        self.global_stats.add_insertion()
        self.stats_by_window[self.get_bin(position)].add_insertion()

    def add_deletion(self, position):
        self.global_stats.add_deletion()
        self.stats_by_window[self.get_bin(position)].add_deletion()

    def get_matches(self):
        return self.global_stats.matches

    def get_mismatches(self):
        return self.global_stats.mismatches
    
    def get_insertions(self):
        return self.global_stats.insertions

    def get_deletions(self):
        return self.global_stats.deletions

    def get_homopolymer_accuracy(self, length):
        return float(self.global_stats.hp_correct[length]) / self.global_stats.hp_count[length]

    def get_identity(self):
        return self.global_stats.get_identity()

    def get_segment_identities(self):
        return [x.get_identity() for x in self.stats_by_window.values()]

    def get_median_segment_identity(self):
        return np.median(self.get_segment_identities())

    def get_bin(self, position):
        return int(position / self.window_size)

def qscore(p):

    if p < 1:
        return -10 * math.log10(1 - p)
    else:
        return 90

def var2key(var):
    return var.CHROM + ":" + str(var.POS) + ":" + var.REF + ">" + str(var.ALT[0])

def key2pos(key):
    fields = key.split(":")
    return int(fields[1])

def load_calls(filename):
    calls = dict()
    vcf_reader = vcf.Reader(open(filename, 'r'))
    for record in vcf_reader:
        key = var2key(record)
        calls[key] = record
    return calls


# convert a pysam record into a pair of aligned strings like so:
#  s1: AGTACTA-GTATAC
#  s2: AGTG-TAGGTATAC
def make_aligned_strings(read, reference_file):

    qseq = read.query_alignment_sequence
    rseq = reference_file.fetch(read.reference_name, read.reference_start, read.reference_end)

    qseq_a_list = list()
    rseq_a_list = list()

    qcurr = 0
    rcurr = 0

    for cigarOp, cigarLength in read.cigar:

        qc = 0
        rc = 0

        # https://www.biostars.org/p/76119/
        if(cigarOp == 0): #match
            qc = cigarLength
            rc = cigarLength
        elif(cigarOp == 1): #insertions
            qc = cigarLength
            rc = 0
        elif(cigarOp == 2): #deletion
            qc = 0
            rc = cigarLength
        elif(cigarOp == 4): #soft clipping
            # ignore since we use query_aligned_sequence
            continue
        elif(cigarOp == 5): # hard clipping
            # ignore
            continue
        else:
            raise ValueError("Unhandled cigar op:" + str(cigarOp))

        # append characters from query/reference to the alignment
        qseq_a_list.append(qseq[qcurr:qcurr+qc])
        rseq_a_list.append(rseq[rcurr:rcurr+rc])

        mc = max(qc, rc)
        qpad = mc - qc
        rpad = mc - rc

        # pad so the lengths match
        qseq_a_list.append("-" * qpad)
        rseq_a_list.append("-" * rpad)

        qcurr += qc
        rcurr += rc

    qaligned = "".join(qseq_a_list).upper()
    raligned = "".join(rseq_a_list).upper()

    return qaligned, raligned

def print_alignment(read, query_aligned, ref_aligned):
    reference_position = read.reference_start
    print("Alignment for %s:%s" % (read.reference_name, reference_position))

    for i in range(0, len(query_aligned), 80):
        print("QUERY    " + query_aligned[i:i+80])
        print("REF      " + ref_aligned[i:i+80] + " " + str(reference_position))
        print("")
        reference_position += len(ref_aligned[i:i+80].replace('-', ''))

def rev_comp_aligned(s):
    trans = str.maketrans("ACGT-", "TGCA-")
    return s.translate(trans)[::-1]

def gather_basic_stats(fp, read, query_aligned, ref_aligned, assembly_accuracy):
    
    reference_position = read.reference_start
    hard_clip_tag = 5
    query_pre_hard_clip = 0
    if read.cigartuples[0][0] == hard_clip_tag:
        query_pre_hard_clip += read.cigartuples[0][1]

    query_post_hard_clip = 0
    if read.cigartuples[-1][0] == hard_clip_tag:
        query_post_hard_clip += read.cigartuples[-1][1]

    # calculate sequence length, including clips    
    sequence_length = query_pre_hard_clip + len(read.seq) + query_post_hard_clip
    
    #
    query_start_position = read.query_alignment_start + query_pre_hard_clip
    query_end_position = read.query_alignment_end + query_pre_hard_clip
    query_position = query_start_position

    # switch strands if rc
    if read.is_reverse:
        query_position = sequence_length - query_end_position
        query_aligned = rev_comp_aligned(query_aligned)
        ref_aligned = rev_comp_aligned(ref_aligned)

    i = 0
    n = len(ref_aligned)
    while i < n:
        if query_aligned[i] == ref_aligned[i]:
            reference_position += 1
            query_position += 1
            i += 1
            assembly_accuracy.add_match(reference_position)

        else:

            # new difference found, find the next matching base
            is_reference_n = False
            j = i
            while query_aligned[j] != ref_aligned[j] and j < n:
                is_reference_n = is_reference_n or ref_aligned[j] == 'N'
                j += 1

            # skip differences at reference gaps
            if is_reference_n:
                i = j
                continue

            # count the error type
            for idx in range(i, j):
                assert(not (query_aligned[idx] == '-' and ref_aligned[idx] == '-'))
                if query_aligned[idx] != '-' and ref_aligned[idx] != '-':
                    assembly_accuracy.add_mismatch(reference_position)

                if ref_aligned[idx] == '-':
                    assembly_accuracy.add_insertion(reference_position)

                if query_aligned[idx] == '-':
                    assembly_accuracy.add_deletion(reference_position)

            # Get difference strings
            q_sub = query_aligned[i:j]
            r_sub = ref_aligned[i:j]

            q_gaps = q_sub.count("-")
            r_gaps = r_sub.count("-")

            offset = 0
            if q_gaps > 0 or r_gaps > 0:
                # append a single base to the start
                q_sub = query_aligned[i-1] + q_sub.replace("-", "")
                r_sub = ref_aligned[i-1] + r_sub.replace("-", "")
                offset = 1
            
            if fp is not None:
                fp.write("%s\t%d\t.\t%s\t%s\t.\tPASS\t.\n" % (read.query_name, query_position - offset + 1, q_sub.upper(), r_sub.upper()))

            # update counters
            reference_position += (j - i - r_gaps)
            query_position += (j - i - q_gaps)
            i = j

def gather_homopolymer_stats(query_aligned, ref_aligned, assembly_accuracy):
    curr_hp_start = 0
    curr_hp_base = 'N'
    n = len(ref_aligned)
    for i in range(0, n):
        if ref_aligned[i] != '-' and ref_aligned[i] != curr_hp_base:
            # end current HP, possibly update stats
            hp_length = i - curr_hp_start
            if hp_length >= args.min_hp_length and hp_length <= args.max_hp_length:
                query_sub = query_aligned[curr_hp_start:i]
                ref_sub = ref_aligned[curr_hp_start:i]
                
                assembly_accuracy.global_stats.hp_count[hp_length] += 1
                assembly_accuracy.global_stats.hp_correct[hp_length] += (query_sub == ref_sub)

            curr_hp_start = i
            curr_hp_base = ref_aligned[i] 
            
parser = argparse.ArgumentParser( description='Calculate the accuracy of a genome assembly by comparing to a reference')
parser.add_argument('--reference', type=str, required=True)
parser.add_argument('--assembly', type=str, required=True)
parser.add_argument('--variants', type=str, required=False)
parser.add_argument('--min-mapping-quality', type=int, required=False, default=0)
parser.add_argument('--min-alignment-length', type=int, required=False, default=50000)
parser.add_argument('--min-hp-length', type=int, required=False, default=4)
parser.add_argument('--max-hp-length', type=int, required=False, default=9)
parser.add_argument('--print-alignment', action='store_true')
parser.add_argument('--print-identity-per-segment', action='store_true')
parser.add_argument('--write-edits', type=str, required=False)
parser.add_argument('--temp-bam', type=str, required=False)
args = parser.parse_args()

out_bam = args.assembly + ".assembly_analysis.sorted.bam"
if args.temp_bam:
    out_bam = args.temp_bam

with open(os.devnull, 'wb') as devnull:
    mm2_cmd = "minimap2 -Y -a -x asm5 %s %s | samtools sort -T %s.tmp -o %s -" % (args.reference, args.assembly, out_bam, out_bam)
    subprocess.check_call(mm2_cmd, stdout=devnull, stderr=devnull, shell=True)

    index_cmd = "samtools index %s" % (out_bam)
    subprocess.check_call(index_cmd, stdout=devnull, stderr=devnull, shell=True)

# Read variants VCF to ignore as errors
variants = dict()
if args.variants is not None:
    variants = load_calls(args.variants)

# Open reference file
reference_file = pysam.FastaFile(args.reference)

# accuracy calculator
assembly_accuracy = AssemblyAccuracy()

edits_fp = None
if args.write_edits is not None:
    edits_fp = open(args.write_edits, "w")

# Calculate the number of matching bases from the alignment
samfile = pysam.AlignmentFile(out_bam)
for read in samfile:

    try:
        if read.mapq < args.min_mapping_quality:
            continue

        if read.query_alignment_length < args.min_alignment_length:
            continue

        query_aligned, ref_aligned = make_aligned_strings(read, reference_file)
        reference_position = int(read.reference_start)

        if args.print_alignment:
            print_alignment(read, query_aligned, ref_aligned)
        
        # accumulate stats for this aligned segment
        gather_basic_stats(edits_fp, read, query_aligned, ref_aligned, assembly_accuracy)
        gather_homopolymer_stats(query_aligned, ref_aligned, assembly_accuracy)

    except ValueError as inst:
        pass

if args.print_identity_per_segment:
    si = assembly_accuracy.get_segment_identities()
    print("index\tqscore")
    for idx, identity in enumerate(si):
        print("%d\t%.2f" % (idx, qscore(identity)))

else:
    identity = assembly_accuracy.get_identity()
    median_identity = assembly_accuracy.get_median_segment_identity()

    # build header and output stats
    header = ["assembly_name", "reference_name", "percent_identity", "qscore", "segment_median_qscore", "num_matches", "num_mismatches", "num_insertions", "num_deletions"]
    stats = [os.path.basename(args.assembly), os.path.basename(args.reference), "%.6f" % (100 * identity), "%.2f" % qscore(identity), "%.2f" % qscore(median_identity),
                                                                                  assembly_accuracy.get_matches(), 
                                                                                  assembly_accuracy.get_mismatches(), 
                                                                                  assembly_accuracy.get_insertions(),
                                                                                  assembly_accuracy.get_deletions()]

    for i in range(args.min_hp_length, args.max_hp_length + 1):
        header.append("%dmer_acc" % (i))
        stats.append("%.3f" % (assembly_accuracy.get_homopolymer_accuracy(i)))

    print("\t".join(header))
    print("\t".join([str(x) for x in stats]))

# clean up temporary files
if not args.temp_bam:
    os.remove(out_bam)
    os.remove(out_bam + ".bai")
