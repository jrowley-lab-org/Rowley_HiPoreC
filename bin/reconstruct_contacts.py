#!/usr/bin/env python3
"""
reconstruct_contacts.py
Read a name-sorted BAM from HiPore-C fragment alignment and reconstruct
multi-way chromatin contacts.

For each original nanopore read (identified by UUID in QNAME):
  - Group all aligned fragment records (QNAME = UUID|fragN|nfTotal)
  - Discard unmapped fragments
  - Apply religation filter: remove the earlier fragment from adjacent
    same-chromosome pairs where gap < 100 bp
    (ported from getpairwise.py lines 34-43)
  - If surviving mapped fragment count >= --min-frags, emit one TSV row
    per fragment

Output TSV columns (no header — header added by EMIT_PARQUET):
  read_id | sample | replicate | n_frags | frag_rank | chrom | start | end | strand | mapq
"""

import pysam
import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Reconstruct HiPore-C multi-way contacts from name-sorted BAM."
    )
    p.add_argument("--bam",       required=True,         help="Name-sorted BAM file")
    p.add_argument("--sample",    required=True,         help="Sample ID (written to output)")
    p.add_argument("--replicate", required=True,         help="Replicate ID (written to output)")
    p.add_argument("--min-frags", type=int, default=2,   help="Min mapped fragments to emit a contact (default: 2)")
    p.add_argument("--out",       required=True,         help="Output TSV file path")
    p.add_argument("--stats",     required=True,         help="Output stats file path")
    return p.parse_args()


def apply_religation_filter(frags):
    """
    Mark and remove the EARLIER fragment from adjacent same-chrom pairs
    where (next.start - prev.end) < 100 bp (likely religations).

    Ported from getpairwise.py: builds indexestopop, then pops in order
    while adjusting for previously removed indices.
    """
    indexes_to_pop = []
    for i in range(len(frags) - 1):
        a, b = frags[i], frags[i + 1]
        if a["chrom"] == b["chrom"] and b["start"] - a["end"] < 100:
            indexes_to_pop.append(i)

    for n_popped, idx in enumerate(indexes_to_pop):
        frags.pop(idx - n_popped)

    return frags


def flush_read(read_id, frags, sample, replicate, min_frags, out_fh, stats):
    if not frags:
        return
    frags.sort(key=lambda f: f["rank"])
    frags = apply_religation_filter(frags)
    n = len(frags)
    stats["n_way"][n] += 1
    if n < min_frags:
        return
    stats["emitted"] += 1
    for new_rank, f in enumerate(frags, start=1):
        out_fh.write(
            f"{read_id}\t{sample}\t{replicate}\t{n}\t{new_rank}\t"
            f"{f['chrom']}\t{f['start']}\t{f['end']}\t{f['strand']}\t{f['mapq']}\n"
        )


def main():
    args = parse_args()
    stats = {
        "total_reads":    0,
        "unmapped_frags": 0,
        "emitted":        0,
        "n_way":          defaultdict(int),
    }

    current_uuid = None
    current_frags = []

    with pysam.AlignmentFile(args.bam, "rb") as bam, open(args.out, "w") as out_fh:
        for read in bam.fetch(until_eof=True):
            parts = read.query_name.split("|")
            uuid  = parts[0]
            rank  = int(parts[1].replace("frag", ""))

            if uuid != current_uuid:
                if current_uuid is not None:
                    flush_read(current_uuid, current_frags, args.sample, args.replicate,
                               args.min_frags, out_fh, stats)
                    stats["total_reads"] += 1
                current_uuid  = uuid
                current_frags = []

            if read.is_unmapped:
                stats["unmapped_frags"] += 1
                continue

            current_frags.append({
                "rank":   rank,
                "chrom":  read.reference_name,
                "start":  read.reference_start,
                "end":    read.reference_end,
                "strand": "-" if read.is_reverse else "+",
                "mapq":   read.mapping_quality,
            })

        if current_uuid is not None:
            flush_read(current_uuid, current_frags, args.sample, args.replicate,
                       args.min_frags, out_fh, stats)
            stats["total_reads"] += 1

    with open(args.stats, "w") as sf:
        sf.write(f"total_reads\t{stats['total_reads']}\n")
        sf.write(f"unmapped_fragments\t{stats['unmapped_frags']}\n")
        sf.write(f"reads_emitted_as_contacts\t{stats['emitted']}\n")
        sf.write("n_way_distribution\n")
        for n in sorted(stats["n_way"]):
            sf.write(f"\t{n}_way\t{stats['n_way'][n]}\n")


if __name__ == "__main__":
    main()
