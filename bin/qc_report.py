#!/usr/bin/env python3
"""
qc_report.py
Aggregate per-sample digest and reconstruction stats files into a
single human-readable QC report.
"""

import argparse
import glob
import os
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate HiPore-C pipeline QC stats.")
    p.add_argument("--stats-dir", required=True, help="Directory containing all *_stats.txt files")
    p.add_argument("--out",       required=True, help="Output report file")
    return p.parse_args()


def parse_stats_file(path):
    """Return dict of key→value from a two-column TSV stats file."""
    data = {}
    section = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("\t"):
                parts = line.strip().split("\t")
                if section and len(parts) == 2:
                    data.setdefault(section, {})[parts[0]] = int(parts[1])
            else:
                parts = line.split("\t")
                if len(parts) == 2:
                    try:
                        data[parts[0]] = int(parts[1])
                    except ValueError:
                        data[parts[0]] = parts[1]
                else:
                    section = parts[0].rstrip(":")
    return data


def main():
    args = parse_args()

    all_files    = glob.glob(os.path.join(args.stats_dir, "*.txt"))
    digest_files = sorted(f for f in all_files if "_digest_stats.txt" in f)
    recon_files  = sorted(f for f in all_files if "_reconstruct_stats.txt" in f)

    with open(args.out, "w") as out:
        out.write("=" * 60 + "\n")
        out.write("HiPore-C Pipeline QC Report\n")
        out.write("=" * 60 + "\n\n")

        out.write("--- DIGESTION STATS ---\n\n")
        for fpath in digest_files:
            sample = os.path.basename(fpath).replace("_digest_stats.txt", "")
            d = parse_stats_file(fpath)
            total      = d.get("total_reads", 0)
            ligation   = d.get("ligation_reads", 0)
            nonlig     = d.get("nonligation_reads", 0)
            lig_pct    = (ligation / total * 100) if total else 0
            out.write(f"Sample: {sample}\n")
            out.write(f"  Total reads       : {total:,}\n")
            out.write(f"  Ligation reads    : {ligation:,}  ({lig_pct:.1f}%)\n")
            out.write(f"  Non-ligation reads: {nonlig:,}\n")
            frag_dist = d.get("fragment_count_distribution", {})
            if frag_dist:
                out.write("  Fragment count distribution:\n")
                for key in sorted(frag_dist, key=lambda k: int(k.split("_")[0])):
                    out.write(f"    {key}: {frag_dist[key]:,}\n")
            out.write("\n")

        out.write("--- CONTACT RECONSTRUCTION STATS ---\n\n")
        for fpath in recon_files:
            sample = os.path.basename(fpath).replace("_reconstruct_stats.txt", "")
            d = parse_stats_file(fpath)
            total    = d.get("total_reads", 0)
            unmapped = d.get("unmapped_fragments", 0)
            emitted  = d.get("reads_emitted_as_contacts", 0)
            emit_pct = (emitted / total * 100) if total else 0
            out.write(f"Sample: {sample}\n")
            out.write(f"  Reads processed           : {total:,}\n")
            out.write(f"  Unmapped fragments skipped: {unmapped:,}\n")
            out.write(f"  Reads emitted as contacts : {emitted:,}  ({emit_pct:.1f}%)\n")
            nway_dist = d.get("n_way_distribution", {})
            if nway_dist:
                out.write("  N-way contact distribution:\n")
                for key in sorted(nway_dist, key=lambda k: int(k.split("_")[0])):
                    out.write(f"    {key}: {nway_dist[key]:,}\n")
            out.write("\n")


if __name__ == "__main__":
    main()
