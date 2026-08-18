nextflow.enable.dsl = 2

// =============================================================================
// CONCAT_FASTQ
// Merge all *.fastq.gz files in a barcode directory into a single file.
// gzip files can be concatenated byte-for-byte and remain valid gzip.
// =============================================================================
process CONCAT_FASTQ {
    tag "${meta.sample_id}_${meta.replicate}"

    input:
    tuple val(meta), path(barcode_dir)

    output:
    tuple val(meta), path("${meta.sample_id}_${meta.replicate}.fastq.gz")

    script:
    """
    cat ${barcode_dir}/*.fastq.gz > ${meta.sample_id}_${meta.replicate}.fastq.gz
    """
}

// =============================================================================
// DIGEST
// Virtual MboI digestion: split reads at GATC sites, tag headers with
// @{UUID}|frag{N}|nf{total} for post-alignment reconstruction.
// Fixes two bugs in the original fastqdigest_mboi.py:
//   - trailing fragment (after last GATC) is now preserved
//   - headers encode parent read ID and fragment rank
// =============================================================================
process DIGEST {
    tag "${meta.sample_id}_${meta.replicate}"
    publishDir { "${params.outdir}/digest/${meta.sample_id}_${meta.replicate}" },
        mode: 'copy', pattern: '*_digest_stats.txt'

    input:
    tuple val(meta), path(fastq_gz)

    output:
    tuple val(meta), path("ligation.fastq.gz"),    emit: ligation
    tuple val(meta), path("nonligation.fastq.gz"), emit: nonligation
    tuple val(meta), path("${meta.sample_id}_${meta.replicate}_digest_stats.txt"), emit: stats

    script:
    """
    digest_mboi.py \\
        -i ${fastq_gz} \\
        -l ligation.fastq.gz \\
        -n nonligation.fastq.gz \\
        -s ${meta.sample_id}_${meta.replicate}_digest_stats.txt
    """
}

// =============================================================================
// ALIGN
// Align digested ligation fragments to the reference with minimap2.
// --secondary=no: suppress secondary alignments (keep one mapping per fragment)
// -F 2048: discard supplementary alignments from chimeric long fragments
// =============================================================================
process ALIGN {
    tag "${meta.sample_id}_${meta.replicate}"

    input:
    tuple val(meta), path(ligation_fastq_gz)

    output:
    tuple val(meta), path("${meta.sample_id}_${meta.replicate}.bam")

    script:
    """
    minimap2 \\
        -ax map-ont \\
        -t ${task.cpus} \\
        --secondary=no \\
        ${params.genome_index} \\
        ${ligation_fastq_gz} \\
    | samtools view -bS -F 2048 \\
        -o ${meta.sample_id}_${meta.replicate}.bam
    """
}

// =============================================================================
// FILTER_MAPQ
// Filter alignments by minimum mapping quality.
// Invoked twice per sample (mapq_low and mapq_high) via channel combine.
// =============================================================================
process FILTER_MAPQ {
    tag "${meta.sample_id}_${meta.replicate}_q${mapq}"

    input:
    tuple val(meta), path(bam), val(mapq)

    output:
    tuple val(meta), path("${meta.sample_id}_${meta.replicate}_q${mapq}.bam"), val(mapq)

    script:
    """
    samtools view -b -q ${mapq} ${bam} \\
        -o ${meta.sample_id}_${meta.replicate}_q${mapq}.bam
    """
}

// =============================================================================
// SORT_INDEX_BAM
// Coordinate-sort and index filtered BAM for genome browser use.
// =============================================================================
process SORT_INDEX_BAM {
    tag "${meta.sample_id}_${meta.replicate}_q${mapq}"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(meta), path(bam), val(mapq)

    output:
    tuple val(meta),
          path("${meta.sample_id}_${meta.replicate}_q${mapq}.sorted.bam"),
          path("${meta.sample_id}_${meta.replicate}_q${mapq}.sorted.bam.bai"),
          val(mapq)

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} \\
        -o ${meta.sample_id}_${meta.replicate}_q${mapq}.sorted.bam
    samtools index ${meta.sample_id}_${meta.replicate}_q${mapq}.sorted.bam
    """
}

// =============================================================================
// NAMESORT_BAM
// Name-sort filtered BAM so that all fragments from the same read are
// contiguous — required by reconstruct_contacts.py for grouping.
// =============================================================================
process NAMESORT_BAM {
    tag "${meta.sample_id}_${meta.replicate}_q${mapq}"

    input:
    tuple val(meta), path(bam), val(mapq)

    output:
    tuple val(meta), path("${meta.sample_id}_${meta.replicate}_q${mapq}.namesorted.bam"), val(mapq)

    script:
    """
    samtools sort -n -@ ${task.cpus} ${bam} \\
        -o ${meta.sample_id}_${meta.replicate}_q${mapq}.namesorted.bam
    """
}

// =============================================================================
// RECONSTRUCT_CONTACTS
// Group fragments by parent UUID, apply religation filter (same-chrom
// adjacent pairs with gap < 100 bp), and emit a per-fragment contact TSV.
// =============================================================================
process RECONSTRUCT_CONTACTS {
    tag "${meta.sample_id}_${meta.replicate}_q${mapq}"

    input:
    tuple val(meta), path(namesorted_bam), val(mapq)

    output:
    tuple val(meta), path("${meta.sample_id}_${meta.replicate}_q${mapq}_contacts.tsv"),
          val(mapq), emit: contacts
    tuple val(meta), path("${meta.sample_id}_${meta.replicate}_q${mapq}_reconstruct_stats.txt"),
          val(mapq), emit: stats

    script:
    """
    reconstruct_contacts.py \\
        --bam ${namesorted_bam} \\
        --sample ${meta.sample_id} \\
        --replicate ${meta.replicate} \\
        --min-frags ${params.min_frags} \\
        --out ${meta.sample_id}_${meta.replicate}_q${mapq}_contacts.tsv \\
        --stats ${meta.sample_id}_${meta.replicate}_q${mapq}_reconstruct_stats.txt
    """
}

// =============================================================================
// EMIT_PARQUET
// Convert the flat per-fragment contact TSV to Apache Parquet (Snappy).
// Schema: read_id | sample | replicate | n_frags | frag_rank |
//         chrom | start | end | strand | mapq
// Query example (DuckDB):
//   SELECT * FROM 'sample_q30.contacts.parquet' WHERE n_frags >= 3
// =============================================================================
process EMIT_PARQUET {
    tag "${meta.sample_id}_${meta.replicate}_q${mapq}"
    publishDir "${params.outdir}/contacts/parquet", mode: 'copy'

    input:
    tuple val(meta), path(contacts_tsv), val(mapq)

    output:
    tuple val(meta),
          path("${meta.sample_id}_${meta.replicate}_q${mapq}.contacts.parquet"),
          val(mapq)

    script:
    """
    python3 - << 'PYEOF'
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

cols   = ['read_id','sample','replicate','n_frags','frag_rank','chrom','start','end','strand','mapq']
dtypes = {'start': 'int32', 'end': 'int32', 'mapq': 'int8', 'n_frags': 'int16', 'frag_rank': 'int16'}
df = pd.read_csv('${contacts_tsv}', sep='\\t', header=None, names=cols, dtype=dtypes)
table = pa.Table.from_pandas(df, preserve_index=False)
pq.write_table(
    table,
    '${meta.sample_id}_${meta.replicate}_q${mapq}.contacts.parquet',
    compression='snappy',
    row_group_size=100_000,
)
PYEOF
    """
}

// =============================================================================
// MERGE_SAMPLE_CONTACTS
// Concatenate all per-replicate contact TSVs for the same sample_id and MAPQ
// threshold into a single merged TSV for sample-level Parquet output.
// =============================================================================
process MERGE_SAMPLE_CONTACTS {
    tag "${sample_meta.sample_id}_q${sample_meta.mapq}"

    input:
    tuple val(sample_meta), path(tsvs)

    output:
    tuple val(sample_meta), path("${sample_meta.sample_id}_q${sample_meta.mapq}_merged.tsv"), val(sample_meta.mapq)

    script:
    """
    cat ${tsvs} > ${sample_meta.sample_id}_q${sample_meta.mapq}_merged.tsv
    """
}

// =============================================================================
// MERGED_EMIT_PARQUET
// Emit a single Parquet file per (sample_id, MAPQ) combining all replicates.
// The replicate column in each row still identifies the source replicate.
// =============================================================================
process MERGED_EMIT_PARQUET {
    tag "${sample_meta.sample_id}_q${sample_meta.mapq}"
    publishDir "${params.outdir}/contacts/parquet/merged", mode: 'copy'

    input:
    tuple val(sample_meta), path(contacts_tsv), val(mapq)

    output:
    tuple val(sample_meta), path("${sample_meta.sample_id}_q${mapq}.merged.contacts.parquet"), val(mapq)

    script:
    """
    python3 - << 'PYEOF'
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

cols   = ['read_id','sample','replicate','n_frags','frag_rank','chrom','start','end','strand','mapq']
dtypes = {'start': 'int32', 'end': 'int32', 'mapq': 'int8', 'n_frags': 'int16', 'frag_rank': 'int16'}
df = pd.read_csv('${contacts_tsv}', sep='\\t', header=None, names=cols, dtype=dtypes)
table = pa.Table.from_pandas(df, preserve_index=False)
pq.write_table(
    table,
    '${sample_meta.sample_id}_q${mapq}.merged.contacts.parquet',
    compression='snappy',
    row_group_size=100_000,
)
PYEOF
    """
}

// =============================================================================
// QC_REPORT
// Aggregate per-sample stats files into one summary report.
// All stats files are mixed into a single channel and staged flat in the work
// dir; qc_report.py detects type by filename suffix to avoid stageAs collisions.
// =============================================================================
process QC_REPORT {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path all_stats

    output:
    path "hiporec_qc_report.txt"

    script:
    """
    qc_report.py --stats-dir . --out hiporec_qc_report.txt
    """
}

// =============================================================================
// Workflow
// =============================================================================
workflow {

    log.info """
        =============================================
        H I P O R E - C   P I P E L I N E  v0.1.0
        =============================================
        samplesheet : ${params.samplesheet}
        outdir      : ${params.outdir}
        genome_index: ${params.genome_index}
        mapq_low    : ${params.mapq_low}
        mapq_high   : ${params.mapq_high}
        min_frags   : ${params.min_frags}
        threads     : ${params.threads}
        """.stripIndent()

    // 1. Parse samplesheet
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .filter { row -> row.sample_id && row.replicate && row.barcode_dir }
        .map { row ->
            def meta = [
                sample_id: row.sample_id,
                replicate: row.replicate,
            ]
            tuple(meta, file(row.barcode_dir, checkIfExists: true))
        }
        .set { samples_ch }

    // 2. Concatenate per-barcode FASTQ files
    CONCAT_FASTQ(samples_ch)

    // 3. Digest at GATC sites
    DIGEST(CONCAT_FASTQ.out)

    // 4. Align ligation fragments
    ALIGN(DIGEST.out.ligation)

    // 5. Filter at two MAPQ thresholds (each BAM × each threshold → two channels)
    ALIGN.out
        .combine(Channel.of(params.mapq_low, params.mapq_high))
        | FILTER_MAPQ

    // 6. Branch: coordinate-sort (for published BAMs) + name-sort (for reconstruction)
    FILTER_MAPQ.out
        .multiMap { meta, bam, mapq ->
            for_sort:     tuple(meta, bam, mapq)
            for_namesort: tuple(meta, bam, mapq)
        }
        .set { filter_split }

    SORT_INDEX_BAM(filter_split.for_sort)
    NAMESORT_BAM(filter_split.for_namesort)

    // 7. Reconstruct multiway contacts from name-sorted BAM
    RECONSTRUCT_CONTACTS(NAMESORT_BAM.out)

    // 8. Emit per-sample Parquet contact files
    EMIT_PARQUET(RECONSTRUCT_CONTACTS.out.contacts)

    // 8b. Merge all replicates per sample_id and emit one Parquet per (sample_id, MAPQ)
    RECONSTRUCT_CONTACTS.out.contacts
        .map { meta, tsv, mapq -> tuple([sample_id: meta.sample_id, mapq: mapq], tsv) }
        .groupTuple()
        | MERGE_SAMPLE_CONTACTS

    MERGED_EMIT_PARQUET(MERGE_SAMPLE_CONTACTS.out)

    // 9. Aggregate QC report — mix all stats into one channel to avoid stageAs collisions
    all_stats_ch = DIGEST.out.stats
        .map { meta, f -> f }
        .mix(RECONSTRUCT_CONTACTS.out.stats.map { meta, f, mapq -> f })
        .collect()

    QC_REPORT(all_stats_ch)

    workflow.onComplete {
        println(workflow.success
            ? "\nDone! Results -> ${params.outdir}\n"
            : "\nPipeline failed. Check error above.\n")
    }
}
