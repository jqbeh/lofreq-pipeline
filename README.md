# lofreq-pipeline

Pipeline to run LoFreq — low-frequency variant detection from Illumina
paired-end reads.

A single Bash script that takes paired FASTQ files and a reference genome,
aligns and deduplicates the reads, prepares the BAM with the quality tags
LoFreq requires, and calls SNVs and indels.

## Usage

```bash
chmod +x lofreq_pipeline.sh
./lofreq_pipeline.sh <accession> <ref.fa> <R1.fastq.gz> <R2.fastq.gz> <threads>
```

Example:

```bash
./lofreq_pipeline.sh SRR12345678 ref/NC_000913.3.fasta \
    reads/SRR12345678_1.fastq.gz reads/SRR12345678_2.fastq.gz 8
```

### Inputs

| Argument | Description |
|---|---|
| `<accession>` | Sample name. Used as the prefix for every output file. |
| `<ref.fa>` | Reference genome FASTA. Indexed by the script with `bwa index`, so its directory must be writable. |
| `<R1.fastq.gz>` | Forward reads. |
| `<R2.fastq.gz>` | Reverse reads. |
| `<threads>` | Threads passed to `bwa mem`. |

### Working directory

The script writes **all output into the directory you run it from**, named
`<accession>_*`. The reference index files (`.amb`, `.ann`, `.bwt`, `.pac`,
`.sa`) are written next to the reference FASTA instead.

Two things follow from this:

- **Run one sample per directory.** Step 6 ends with `rm -f *.sam`, which
  removes every SAM file in the working directory, not only this sample's.
- **Don't run two samples against the same reference at the same time.** Each
  run calls `bwa index` on the reference, so concurrent runs overwrite the same
  index files while the other is reading them.

### Requirements

`bwa`, `samtools`, `lofreq` and `bamtools` must be on `PATH`.

The script activates a Conda environment by absolute path:

```bash
source /home/jbbeh/miniconda3/etc/profile.d/conda.sh
conda activate everything
```

**Edit these two lines** to point at your own Conda installation and
environment, or delete them and activate the environment yourself before
running.

## Workflow

### 1. Alignment

| Step | Command | Purpose |
|---|---|---|
| 1 | `bwa index` | Index the reference genome. |
| 2 | `bwa mem` | Align the paired reads, producing a SAM file. |

### 2. Post-processing

| Step | Command | Purpose |
|---|---|---|
| 3 | `samtools sort -n` | Sort by read name so mates sit together. |
| 4 | `samtools fixmate -m` | Add mate-score (`ms`) tags used by `markdup`. |
| 5 | `samtools sort` | Sort by genomic coordinate. |
| 6 | `samtools markdup -r -S` | Remove PCR/optical duplicates, then delete the intermediate SAM files. |
| 7 | `samtools flagstat` | Summary mapping statistics. |
| 8 | `samtools depth` | Per-base coverage. |

Duplicate removal is optional but recommended for Illumina paired-end data at
roughly 50–200× depth. Duplicates are reads originating from a single DNA
fragment; the equivalent tool in GATK is Picard `MarkDuplicates`.

### 3. Variant calling

| Step | Command | Purpose |
|---|---|---|
| 9 | `bamtools index` | Index the deduplicated BAM. |
| 10 | *(check)* | Report whether `BI`/`BD` indel-quality tags are already present. Diagnostic only — step 12 runs either way. |
| 11 | `lofreq viterbi` | Probabilistic realignment to correct mapping errors, then re-sort. |
| 12 | `lofreq indelqual --dindel` | Insert indel qualities (`BI` = insertion, `BD` = deletion). Usually needed for bacterial and viral paired-end data. |
| 13 | `lofreq alnqual -b` | Add base alignment qualities (BAQ), including for indels. |
| 14 | `lofreq call --call-indels` | Call SNVs and indels. |

Step 13 matters: `lofreq call` adds BAQ on its own but **not** indel alignment
qualities, so `alnqual` has to run first.

Two options are left commented out at step 14:

- `--no-default-filter` — keep all calls and filter them yourself downstream.
- `-a` — p-value threshold (default `0.01`, i.e. Q20).

## Outputs

All prefixed with `<accession>`:

| File | Contents |
|---|---|
| `_aln_fixmate.bam` | After `fixmate`. |
| `_aln_sorted.bam` | Coordinate-sorted, before deduplication. |
| `_aln_sorted_dedup.bam` (`.bai`) | Deduplicated alignment. |
| `_mapping_stats.txt` | `samtools flagstat` output. |
| `_mapping_depth.txt` | Per-base depth. |
| `_lofreq_realign.bam` (`.bai`) | After `lofreq viterbi`. |
| `_lofreq.indelqual.bam` (`.bai`) | After `lofreq indelqual`. |
| `_lofreq_realign_alnq.bam` (`.bai`) | After `lofreq alnqual` — the BAM used for calling. |
| `_lofreq.vcf` | **Final variant calls.** |

The intermediate SAM files are deleted at step 6. The intermediate BAMs are
kept; delete them once you have the VCF if disk space matters.

## Note on read filtering

Do not pre-filter reads aggressively on base or mapping quality. LoFreq's model
is built to handle quality variation, and filtering introduces the calling bias
it was designed to be robust against.
