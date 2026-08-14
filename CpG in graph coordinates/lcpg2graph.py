

import os
import sys
import json
import copy
import time
import argparse
import subprocess
import multiprocessing







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



def parse_gfa(ggfp1, ggfp2, selected_chrom, sample_id, haplo_id):
    segments = {}
    all_walks = {}
    with open(ggfp2) as ggfh:
        for l in ggfh:
            if not l.startswith("W"):
                continue

            l = l.strip().split("\t")

            if l[1] != sample_id:
                continue

            if l[2] != haplo_id:
                continue

            chrom = l[3]
            start = int(l[4])
            end = int(l[5])

            if chrom != selected_chrom:
                continue

            walk_str = l[6]
            walk_list = walk2list(walk_str)

            all_walks[(chrom, start, end)] = walk_list
            for ori, seg_id in walk_list:
                segments[seg_id] = None
            # print(chrom, start, end)
            # break

    with open(ggfp1) as ggfh:
        for l in ggfh:
            if not l.startswith("S"):
                continue

            l = l.strip().split("\t")

            seg_id = l[1]
            seq = l[2]

            if seg_id not in segments:
                continue
            segments[seg_id] = seq

    for seg_id, seq in segments.items():
        assert seq != None

    return segments, all_walks






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
            assert res >=1 and res <=2
            end = start + 2

            if chrom not in peaks:
                peaks[chrom] = []

            peaks[chrom].append((start, end, tuple(anno)))

    # Need to sort peaks by start position for future use
    for chrom in peaks:
        peaks[chrom].sort(key=lambda x: x[0])

    return peaks


def peaks_all_chr(peak_table_fp):
    res = set()
    with open(peak_table_fp) as peak_table_fh:
        for l in peak_table_fh:
            l = l.strip().split("\t")
            chrom = l[0]
            res.add(chrom)
    return res


def sequence_reverse_complement(seq):
    return seq.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def linear_range_to_graph_coord(segments, all_walks, peaks, selected_chrom, debug_flag=True):
    peak_selected = peaks[selected_chrom]

    start_segment_index = 0
    peak_index = 0

    intended_walk = None
    for chrom, start, end in all_walks.keys():
        if chrom != selected_chrom:
            continue
        intended_walk = all_walks[(chrom, start, end)]
    assert intended_walk != None

    segment_range_start = 0
    segment_range_end = len(segments[intended_walk[0][1]])

    result_lines = []

    while peak_index < len(peak_selected):

        peak_start, peak_end, peak_annotation = peak_selected[peak_index]
        peak_width = peak_end - peak_start

        if not segment_range_start <= peak_start < segment_range_end:
            start_segment_index += 1
            segment_range_start = segment_range_end
            segment_range_end += len(segments[intended_walk[start_segment_index][1]])
            continue

        # print(peak_index, start_segment_index, peak_start)

        peak_path_indexes = [start_segment_index]
        peak_this_segment_index = start_segment_index
        peak_this_segment_start = peak_start - segment_range_start
        peak_rest = peak_width - (len(segments[intended_walk[start_segment_index][1]]) - peak_this_segment_start)

        # I know this code looks stupid... For a very big segment which completely includes the peak, it needs to be handled separately
        if peak_rest > 0:
            while True:
                # print(peak_this_segment_index, peak_rest)
                peak_this_segment_index += 1
                peak_rest_tmp = peak_rest - len(segments[intended_walk[peak_this_segment_index][1]])
                peak_path_indexes.append(peak_this_segment_index)

                if peak_rest_tmp <= 0:
                    break
                peak_rest = peak_rest_tmp

            peak_this_segment_end = peak_rest
        else:
            peak_this_segment_end = peak_this_segment_start + peak_width

        peak_graph_path_last_segment_sequence_length = len(segments[intended_walk[peak_path_indexes[-1]][1]])
        assert 0 < peak_this_segment_end <= peak_graph_path_last_segment_sequence_length

        peak_path = []
        peak_path_len = 0
        peak_path_str = ""
        peak_sequence = ""

        for i in peak_path_indexes:
            peak_path.append(intended_walk[i])
            peak_path_str += {"+": ">", "-": "<"}[intended_walk[i][0]] + intended_walk[i][1]

            if intended_walk[i][0] == "+":
                peak_sequence += segments[intended_walk[i][1]]
            else:
                peak_sequence += sequence_reverse_complement(segments[intended_walk[i][1]])
            peak_path_len += len(segments[intended_walk[i][1]])

        peak_sequence = peak_sequence[peak_this_segment_start:peak_this_segment_start + peak_width]

        # Recalculate peak width on graph, make sure it is correct
        pl = -1
        for i in peak_path_indexes:
            seg = intended_walk[i][1]
            seq = segments[seg]
            seg_length = len(seq)

            seg_start = 0
            seg_end = seg_length - 1

            if i == peak_path_indexes[0]:
                seg_start = peak_this_segment_start
            if i == peak_path_indexes[-1]:
                seg_end = peak_this_segment_end
            pl += seg_end - seg_start + 1

        pld = abs(pl - peak_width)
        assert pld == 0

        res_line = [selected_chrom, peak_start, peak_end] + list(peak_annotation) +[peak_path_str, peak_this_segment_start,
                    peak_this_segment_end, peak_path_len]
        debug_info = [peak_sequence]

        if debug_flag:
            res_line += debug_info

        rls = "\t".join(map(str, res_line))
        result_lines.append(rls)

        debugx = [peak_width]

        # print(peak_start, peak_end, peak_annotation)
        # print(segment_range_start, segment_range_end)
        # print(rls)
        # print(debug)
        # print()

        peak_index += 1
        # break

        # if len(result_lines) > 100:
        #    break

        # if pld >= 1:
        #    print(res_line[:-1], pl)
        # print(pld)

    return result_lines


def mpfunc(ggfp1, ggfp2, peak_table_fp, selected_chrom, sample_id, haplo_id, debug_flag):
    print(f"Processing {selected_chrom}...", file=sys.stderr)
    segments, all_walks = parse_gfa(ggfp1, ggfp2, selected_chrom, sample_id, haplo_id)
    peaks = read_bed(peak_table_fp, selected_chrom=selected_chrom)
    print(f"Processing {len(peaks[selected_chrom])} CpGs...", file=sys.stderr)
    result_lines = linear_range_to_graph_coord(segments, all_walks, peaks, selected_chrom, debug_flag=debug_flag)
    print(f"Finished processing.", file=sys.stderr)

    return result_lines



if __name__ == "__main__":

    threads = 10
    debug = True

    args = sys.argv[1:]
    if len(args) >= 4:
        ggfp1, ggfp2, cpg_table_fp, sample_id, haplo_id = args[:5]
        if len(args) >= 6:
            threads = int(args[5])
        if len(args) >= 7:
            debug = False

    else:
        # print("Usage: python lr2gc.py path_to_genome_graph path_to_bed selected_chrom sample_id haplo_id [threads] [debug] > output.txt")
        sys.exit(1)


    """
    ggfp = f"/scratch/wzhang/ref/graph/mouse/all-crosses/all-crosses.full.gfa"
    peak_table_fp = f"./CC001_ids.txt"
    
    selected_chrom = "chr1"
    sample_id = "CC001"
    haplo_id = "0"
    
    debug = True
    """
    print(ggfp1, ggfp2, cpg_table_fp, sample_id, haplo_id, threads, debug, file=sys.stderr)

    all_chr = peaks_all_chr(cpg_table_fp)

    pool = multiprocessing.Pool(threads)
    results = []
    for selected_chrom in all_chr:
        results.append(
            pool.apply_async(mpfunc, (ggfp1, ggfp2, cpg_table_fp, selected_chrom, sample_id, haplo_id, debug)))

    pool.close()
    pool.join()
    print(f"Pool closed and joined.", file=sys.stderr)

    # print(f"Writing results to {out_fp}...")

    for r in results:
        try:
            for l in r.get():
                print(l)
                # break
        except:
            continue







