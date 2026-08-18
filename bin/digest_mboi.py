#!/usr/bin/env python3
"""
digest_mboi.py
Virtual digestion of HiPore-C long nanopore reads at MboI (GATC) sites.

Fixes over the original fastqdigest_mboi.py:
  1. Header tagging: each fragment's QNAME encodes parent UUID, rank, and
     total fragment count as @{UUID}|frag{N}|nf{total}.  Splitting the QNAME
     on '|' after alignment recovers the parent read, enabling multi-way
     contact reconstruction.
  2. Trailing fragment preserved: the original script silently dropped the
     last genomic fragment of every multi-site read.  This version collects
     all fragments before writing so the trailing piece is included.
  3. Gzip output for space efficiency (minimap2 reads .fastq.gz natively).
  4. Stats file: fragment-count distribution written to --stats.

Cutting convention (preserved from original fastqdigest_mboi.py):
  On finding GATC at position p in the read sequence:
    - Current fragment  = seq[current_pos : p+1]  (includes 'G')
    - Discarded region  = seq[p+1 : p+5]          ('ATC' + 1 base)
    - Next fragment starts at p+5
"""

import gzip
import argparse
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Virtual MboI digestion of HiPore-C FASTQ with fragment tracking.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("-i", "--fastqfile",      dest="fastq_file",      required=True, metavar="FILE",
                   help="Input .fastq.gz file")
    p.add_argument("-l", "--ligationout",    dest="ligation_out",     required=True, metavar="FILE",
                   help="Output .fastq.gz: reads with GATC junctions (one entry per fragment)")
    p.add_argument("-n", "--nonligationout", dest="nonligation_out",  required=True, metavar="FILE",
                   help="Output .fastq.gz: reads with no GATC junction (one entry per read)")
    p.add_argument("-s", "--stats",          dest="stats_out",        required=True, metavar="FILE",
                   help="Output stats text file")
    return p.parse_args()


def split_at_gatc(seq, qual):
    """
    Return list of (frag_seq, frag_qual) or None if no GATC found.

    Matches the cutting convention of the original fastqdigest_mboi.py:
      GATC at position p → fragment ends at p+1 (includes 'G'),
      next fragment starts at p+5 (skips 'ATC' + one subsequent base).
    """
    fragments = []
    current_pos = 0
    search_from = 0
    n = len(seq)

    while True:
        pos = seq.find("GATC", search_from)
        if pos == -1 or pos > n - 4:
            break
        cut = pos + 1
        fragments.append((seq[current_pos:cut], qual[current_pos:cut]))
        current_pos = pos + 5
        search_from = pos + 5

    if not fragments:
        return None

    # trailing fragment after the last GATC site
    if current_pos < n:
        fragments.append((seq[current_pos:], qual[current_pos:]))

    return fragments


def iter_fastq_gz(path):
    """Yield (header_line, seq, qual) from a gzipped FASTQ file."""
    with gzip.open(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().rstrip("\n")
            fh.readline()               # '+' separator
            qual = fh.readline().rstrip("\n")
            yield header.rstrip("\n"), seq, qual


def extract_uuid(header_line):
    """Return the nanopore read UUID (first whitespace-delimited token after '@')."""
    return header_line.split()[0].lstrip("@")


def main():
    args = parse_args()
    frag_dist = defaultdict(int)
    total_reads = ligation_reads = nonligation_reads = 0

    with gzip.open(args.ligation_out, "wt") as lig_fh, \
         gzip.open(args.nonligation_out, "wt") as nonlig_fh:

        for header, seq, qual in iter_fastq_gz(args.fastq_file):
            total_reads += 1
            uuid = extract_uuid(header)
            fragments = split_at_gatc(seq, qual)

            if fragments is None:
                nonligation_reads += 1
                frag_dist[1] += 1
                nonlig_fh.write(f"@{uuid}|frag1|nf1\n{seq}\n+\n{qual}\n")
            else:
                ligation_reads += 1
                n = len(fragments)
                frag_dist[n] += 1
                for rank, (fseq, fqual) in enumerate(fragments, start=1):
                    lig_fh.write(f"@{uuid}|frag{rank}|nf{n}\n{fseq}\n+\n{fqual}\n")

    with open(args.stats_out, "w") as sf:
        sf.write(f"total_reads\t{total_reads}\n")
        sf.write(f"ligation_reads\t{ligation_reads}\n")
        sf.write(f"nonligation_reads\t{nonligation_reads}\n")
        sf.write("fragment_count_distribution\n")
        for n_frags in sorted(frag_dist):
            sf.write(f"\t{n_frags}_frags\t{frag_dist[n_frags]}\n")


if __name__ == "__main__":
    main()
