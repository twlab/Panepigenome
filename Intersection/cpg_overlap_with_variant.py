
import os
import sys
import json
import copy
import time
import argparse
import subprocess
import multiprocessing







def sequence_reverse_complement(seq):
    return seq.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def walk2list(walk_str):
    walk = []

    ori = None
    seg_id = ""


    for c in walk_str:

        if c in "<>":

            if ori != None:
                walk.append((ori, seg_id))
                ori = ""
                seg_id = ""

            if c == ">":
                ori = "+"
            else:
                ori = "-"
            continue

        seg_id += c

    walk.append((ori, seg_id))

    reconstructed_walk_str = ""
    for ori, seg_id in walk:
        reconstructed_walk_str += {"+":">", "-":"<"}[ori] + seg_id
    assert reconstructed_walk_str == walk_str


    return walk



def parse_gfa(ggfp1):
    segments = {}

    with open(ggfp1) as ggfh:
        for l in ggfh:
            if not l.startswith("S"):
                continue

            l = l.strip().split("\t")

            seg_id = l[1]
            seq = l[2]

            segments[seg_id] = len(seq)

    return segments






def read_bed(peak_table_fp, selected_chrom=None):
    peaks = {}
    with open(peak_table_fp) as peak_table_fh:

        for l in peak_table_fh:
            l = l.strip().split("\t")
            chrom, start, end, *anno = l
            start = int(start)
            end = int(end)

            if selected_chrom is not None and chrom != selected_chrom:
                continue

            if end < start:
                continue

            res = end - start
            # assert res >=1 and res <=2
            # end = start + 2

            if chrom not in peaks:
                peaks[chrom] = []

            peaks[chrom].append((start, end, tuple(anno)))

    # Need to sort peaks by start position for future use
    for chrom in peaks:
        peaks[chrom].sort(key=lambda x: x[0])

    return peaks





def overlap_calc(w1, w2):
    overlap = 0

    s1 = {}
    s2 = {}
    for ori, seg_id, ss, se in w1:
        s1[seg_id] = (ss, se)
    for ori, seg_id, ss, se in w2:
        s2[seg_id] = (ss, se)

    overlap_segments = set(s1.keys()).intersection(set(s2.keys()))

    if len(overlap_segments) == 0:
        return overlap


    for seg_id in overlap_segments:
        ss1, se1 = s1[seg_id]
        ss2, se2 = s2[seg_id]

        o = min([se1, se2]) - max([ss1, ss2])
        if o <= 0:
            continue
        overlap += o

    return overlap



def variants_reader(fp):
    res = []
    with open(fp) as f:
        for l in f:
            l = l.strip().split("\t")
            if len(l) == 19:
                l.append("-")

            assert len(l) == 20
            res.append(l)
            # break

    return res


def cpg_reader(fp):
    res = []
    with open(fp) as f:
        for l in f:
            l = l.strip().split("\t")
            assert len(l) == 5
            res.append(l)
            # break

    return res








argv = sys.argv[1:]
assert len(argv) == 3

graph_fp, f1, f2 = argv
# print(graph_fp)
# print(f1)
# print(f2)


segments = parse_gfa(graph_fp)


variant_info = variants_reader(f1)
cpg_info = cpg_reader(f2)

variants_by_segments = {}

for vi,v in enumerate(variant_info):
    path_str = v[15]
    path = walk2list(path_str)

    for ori, seg_id in path:
        if seg_id not in variants_by_segments:
            variants_by_segments[seg_id] = []
        variants_by_segments[seg_id].append(vi)



def graph_range_to_bp(graph_range):
    res = set()


    path = walk2list(graph_range[0])
    first_seg_start = int(graph_range[1])
    last_seg_end = int(graph_range[2])

    for si, (ori, seg_id) in enumerate(path):
        seg_start = 0
        seg_end = segments[seg_id]-1

        if si == 0:
            seg_start = first_seg_start
        if si == len(path) - 1:
            seg_end = last_seg_end-1

        if ori == "-":
            seg_start, seg_end = segments[seg_id] - seg_end - 1, segments[seg_id] - seg_start - 1
        # print(seg_id, seg_start, seg_end)

        for bp in range(seg_start, seg_end+1):
            res.add((seg_id, bp))

    return res

def cpg_overlap_with_variant(cpg_loc, variant_loc):
    overlap = 0

    cpg_bp_locations = graph_range_to_bp(cpg_loc)
    variant_bp_locations = graph_range_to_bp(variant_loc)
    # print(cpg_loc, cpg_bp_locations)
    # print()
    overlap = len(cpg_bp_locations.intersection(variant_bp_locations))

    return overlap


result = set()
for cpg_line_num, cpg in enumerate(cpg_info):

    path = walk2list(cpg[0])
    interesting_variants = []

    for ori, seg_id in path:
        interesting_variants += variants_by_segments.get(seg_id, [])

    interesting_variants = list(sorted(set(interesting_variants)))

    if len(interesting_variants) == 0:
        continue


    for ivln in interesting_variants:
        cpg_loc = cpg[0:3]
        variant_loc = variant_info[ivln][15:18]
        ol = cpg_overlap_with_variant(cpg_loc, variant_loc)
        if ol > 0:
            result.add((cpg_line_num, ivln, ol))

result = list(result)


for r in result:
    cpg_line_num, ivln, ol = r

    l = [cpg_line_num, ivln, ol] + cpg_info[cpg_line_num] + variant_info[ivln]
    l = [str(x) for x in l]
    print("\t".join(l))
















