# HiPore-C Nextflow Pipeline Plan

## Context
Build a containerized Nextflow (DSL2) pipeline in `/Zulu/timrez/Pipelines/Rowley_HiPoreC/` to process HiPore-C nanopore data. The existing ad-hoc workflow at `/Tango/RawNanopore/HiPoreCAnalysis/HiPoreC_pipeline.bash` is not reproducible, cannot handle multiple samples/replicates in one run, uses Juicer `.hic` as its only output, and has two bugs in the digestion script. This pipeline replaces it with a formal multi-sample workflow producing multiway contact files at two MAPQ thresholds (≥1 and ≥30), with no .hic or .mcool output.

**Key discoveries from existing code:**
- Reference genome: `/home/genomefiles/human/hg38/hg38.fa`, pre-built minimap2 index at `/home/genomefiles/human/hg38/minimap2/hg38.mmi`
- Containerization: Singularity (not Docker); Nextflow 24.10.4 installed at `/Zulu/timrez/Pipelines/nextflow`
- Existing DSL2 Nextflow pattern to follow: `/Zulu/timrez/Pipelines/chipseq-nf/main.nf`
- Religation filter logic to port from: `/Tango/RawNanopore/HiPoreCAnalysis/getpairwise.py`

---

## Bugs to Fix in Digestion Script

The `bin/digest_mboi.py` (adapted from `fastqdigest_mboi.py`) must fix two bugs:

1. **Last fragment lost**: The sliding-window loop exits without writing the trailing fragment after the final GATC. Fix: collect all `(seq, qual)` tuples per read in a list first, then write all at once.

2. **No read tracking in headers**: All fragments of one read get an identical header, making it impossible to group them after alignment. Fix: rewrite each fragment header as `@{UUID}|frag{N}|nf{total}` where UUID is the original nanopore read ID (first whitespace-delimited token after `@`), N is the 1-based fragment rank, and total is the count of fragments from that read.

Also add: `-s/--stats` argument to write a stats file (total reads, ligation reads, fragment count distribution).

---

## Directory Structure

```
/Zulu/timrez/Pipelines/Rowley_HiPoreC/
├── main.nf                      # Single-file DSL2 pipeline (all processes inline)
├── nextflow.config              # Profiles, params, Singularity config
├── Samplesheets/
│   └── example_samplesheet.csv
├── bin/                         # Auto-added to PATH by Nextflow
│   ├── digest_mboi.py           # Fixed digestion script
│   ├── reconstruct_contacts.py  # BAM → per-fragment contact TSV
│   ├── tsv_to_jsonl.py          # TSV → JSON-Lines walks (if Format 3 chosen)
│   └── qc_report.py             # Collects stats → summary report
└── containers/
    └── hiporec.def              # Singularity definition
```

---

## Samplesheet Format

```csv
sample_id,replicate,timepoint,barcode_dir
SFB_T0,rep1,T0,/Tango/RawNanopore/backup_Dec_17_2024/Puspo/HiPoreC_01_29_2024/barcode01
SFB_T0,rep2,T0,/Tango/RawNanopore/backup_Dec_17_2024/Puspo/HiPoreC_N_01_29_2024/barcode01
SFB_T25,rep1,T25,/Tango/RawNanopore/backup_Dec_17_2024/Puspo/HiPoreC_01_29_2024/barcode02
```

`barcode_dir` points to the directory; all `*.fastq.gz` files within it are concatenated automatically.

---

## Pipeline DAG

```
samplesheet (CSV)
    ↓ CONCAT_FASTQ    — cat all *.fastq.gz per barcode dir
    ↓ DIGEST          — split at GATC, tag headers @UUID|fragN|nfTotal
    ↓ ALIGN           — minimap2 -ax map-ont --secondary=no → BAM
    ↓ FILTER_MAPQ     — samtools view -q {1,30} (two parallel branches)
    ├─ SORT_INDEX_BAM — coordinate sort + index → published BAM
    └─ NAMESORT_BAM   — name sort → input for reconstruction
    ↓ RECONSTRUCT     — group by UUID, filter religations, emit contact TSV
    ↓ EMIT_FORMAT     — convert to chosen output format(s)
    ↓ QC_REPORT       — collect all stats → summary
```

---

## Key Parameters (`nextflow.config`)

```groovy
params {
    samplesheet    = "${projectDir}/Samplesheets/example_samplesheet.csv"
    outdir         = "${projectDir}/results"
    genome_index   = "/home/genomefiles/human/hg38/minimap2/hg38.mmi"
    threads        = 8
    mapq_low       = 1
    mapq_high      = 30
    output_format  = "parquet"   // chosen: Parquet only
    min_frags      = 2           // minimum mapped fragments to emit a contact
}

process.container   = '/Zulu/timrez/Containers/hiporec.sif'
singularity.enabled = true
singularity.autoMounts = true
```

---

## Process Summaries

| Process | Input | Output | Tool |
|---|---|---|---|
| CONCAT_FASTQ | barcode_dir | merged.fastq.gz | cat |
| DIGEST | merged.fastq.gz | ligation.fastq, nonligation.fastq, stats | digest_mboi.py |
| ALIGN | ligation.fastq | raw.bam | minimap2 |
| FILTER_MAPQ | raw.bam × {1,30} | filtered_q{N}.bam | samtools view -q |
| SORT_INDEX_BAM | filtered.bam | sorted.bam + .bai | samtools sort/index |
| NAMESORT_BAM | filtered.bam | namesorted.bam | samtools sort -n |
| RECONSTRUCT_CONTACTS | namesorted.bam | contacts.tsv, stats | reconstruct_contacts.py |
| EMIT_FCT / EMIT_PARQUET / EMIT_JSONL | contacts.tsv | output format | bgzip+tabix / pandas+pyarrow / python json |
| QC_REPORT | all stats files | qc_report.txt | qc_report.py |

**RECONSTRUCT_CONTACTS logic** (`bin/reconstruct_contacts.py`):
- Iterate name-sorted BAM; group records by UUID (split QNAME on `|`)
- Skip unmapped fragments; port religation filter from `getpairwise.py` (adjacent same-chrom fragments < 100 bp gap → collapse/skip)
- If `len(mapped_fragments) >= min_frags`: emit one TSV row per fragment: `read_id | sample | replicate | n_frags | frag_rank | chrom | start | end | strand | mapq`

---

## Three Output Format Options

### Format 1: Fragment-indexed BGZ+tabix table (`.fct.tsv.bgz` + `.tbi`)
One row per fragment. Sorted by chrom+start, bgzip-compressed, tabix-indexed.
- **Genomic lookup**: `tabix sample_q30.fct.tsv.bgz chr1:10000-20000` — fastest random-access genomic queries
- **Directionality**: `frag_rank` column preserves fragment order along the original read; reconstruct full walk by `GROUP BY read_id ORDER BY frag_rank`
- **Filter by n contacts**: `awk '$4 >= 3'` on the `n_frags` column before tabix, or filter post-query
- **Tooling**: tabix, bedtools, pysam — standard bioinformatics CLI

### Format 2: Apache Parquet columnar table (`.contacts.parquet`)
Identical schema to Format 1, columnar binary format with Snappy compression.
- **Genomic lookup**: Partition pruning by chrom + predicate pushdown — fast for `WHERE chrom='chr1' AND start BETWEEN x AND y`
- **Directionality**: `frag_rank` column; pandas/polars sort is native and fast
- **Filter by n contacts**: `df[df.n_frags >= 3]` or `SELECT * FROM '*.parquet' WHERE n_frags >= 3` (DuckDB, no data load)
- **Tooling**: pandas, polars, DuckDB, R/arrow — best for interactive analysis notebooks

### Format 3: JSON-Lines read walks (`.walks.jsonl.gz`)
One JSON object per read, fragments stored as an ordered array.
```json
{"read_id":"82c21677...","sample":"SFB_T0","replicate":"rep1","n_frags":4,
 "fragments":[{"rank":1,"chrom":"chr1","start":13718,"end":14200,"strand":"+","mapq":45},
              {"rank":2,"chrom":"chr1","start":200702,"end":201100,"strand":"-","mapq":42},
              {"rank":3,"chrom":"chr3","start":56001234,"end":56002100,"strand":"+","mapq":38},
              {"rank":4,"chrom":"chr7","start":12345678,"end":12346200,"strand":"-","mapq":30}]}
```
- **Genomic lookup**: No index — requires streaming scan; not suitable for ad-hoc range queries
- **Directionality**: Full ordered walk is the atomic unit; no reconstruction needed — array index IS the fragment order on the DNA molecule
- **Filter by n contacts**: `jq 'select(.n_frags >= 3)'` on stream — trivial
- **Tooling**: Python json/ijson, jq, R jsonlite — best for higher-order/graph analyses where the full N-way hyperedge is the input

**Combinations are possible**: Format 1 or 2 for lookup + Format 3 for graph analysis is a natural pairing. Format 1 + Format 2 is redundant (same data, different access patterns).

---

## Output Directory Layout

```
results/
├── digest/{sample_replicate}/ligation.fastq, nonligation.fastq, digest_stats.txt
├── bam/{sample_replicate}_q{1,30}.sorted.bam[.bai]
├── contacts/
│   ├── fct/{sample_replicate}_q{1,30}.fct.tsv.bgz[.tbi]
│   ├── parquet/{sample_replicate}_q{1,30}.contacts.parquet
│   └── jsonl/{sample_replicate}_q{1,30}.walks.jsonl.gz
└── qc/hiporec_qc_report.txt
```

---

## Singularity Container

`containers/hiporec.def` — Ubuntu 22.04 base with: minimap2 2.28, samtools 1.21, htslib 1.21 (bgzip+tabix), Python 3 + pysam + pandas + pyarrow.

Build: `sudo singularity build /Zulu/timrez/Containers/hiporec.sif containers/hiporec.def`

---

## Verification

1. **Algorithm correctness**: Run `bin/digest_mboi.py` on one small `.fastq.gz`, verify headers have `|fragN|nfTotal`, trailing fragment is present, fragment base counts sum correctly.
2. **Single-sample test**: Run pipeline on barcode01 only; check `samtools flagstat` on BAM, `tabix` query on FCT, `jq` on JSONL first line.
3. **Reconstruction sanity**: Count n-way distribution from output — expect decreasing counts at 2→3→4→5+ fragments per read.
4. **Regression**: Compare pairwise contact count against existing pipeline output for same barcode.
