#!/usr/bin/env bash
set -euo pipefail

########################################
# Usage:
# chmod +x lofreq_pipeline.sh
# ./lofreq_pipeline.sh <accession> <ref.fa> <R1.fastq.gz> <R2.fastq.gz> <threads>

# All output files will be written in the working dir.
########################################


# Define input
SAMPLE="$1"
REF="$2"
R1="$3"
R2="$4"
THREADS="$5"

source /home/jbbeh/miniconda3/etc/profile.d/conda.sh

set +u
conda activate everything
set -u

# 1. Align PE reads to ref
# index your ref
printf '\n[Step 1] Indexing reference genome...\n'
bwa index ${REF}
# align reads to ref
printf '\n[Step 2] Aligning reads...\n'
bwa mem -t ${THREADS} ${REF} ${R1} ${R2} > ./${SAMPLE}_aln.sam

# 2. Post-processing for bwa-mem
# sort alignment
# this makes it easier for downstream analysis 
# to find mate pairs
printf '\n[Step 3] Sorting alignment...\n'
samtools sort -n -O sam ${SAMPLE}_aln.sam > ${SAMPLE}_aln_sorted.sam
  # -n: sort by read name rather than genomic coordinate
  # -O sam: output in sam format
# add mate information
printf '\n[Step 4] Adding mate information...\n'
samtools fixmate -m -O bam ${SAMPLE}_aln_sorted.sam ${SAMPLE}_aln_fixmate.bam
  # -m: Add ms (mate score) tags. These are used by markdup (below) to select the best reads to keep.
  # -O bam: outputs bam format
# sort bam
printf '\n[Step 5] Sorting bam file...\n'
samtools sort -O bam ${SAMPLE}_aln_fixmate.bam > ${SAMPLE}_aln_sorted.bam
  # -O bam: output in bam format
# remove duplicate reads
# this step is optional but recommended for snp calling
# on ILL PE reads with depth >50x and <200x.
# Duplicate reads are defined as originating from a single fragment of DNA.
# Similar tool is MarkDuplicates (Picard) - GATK
printf '\n[Step 6] Removing duplicate reads...\n'
samtools markdup -r -S ${SAMPLE}_aln_sorted.bam ${SAMPLE}_aln_sorted_dedup.bam
  # Markdup can simply *mark* the duplicate reads
  # -r: tells it to remove those reads
  # -S: mark supplementary mappings, but does not remove
# remove sam files to save disk space
rm -f *.sam

# Get mapping statistics
printf '\n[Step 7] Get mapping statistics...\n'
samtools flagstat ${SAMPLE}_aln_sorted_dedup.bam > \
    ${SAMPLE}_mapping_stats.txt
# get more detailed stats on bam coverage
printf '\n[Step 8] Get mapping coverage...\n'
samtools depth ${SAMPLE}_aln_sorted_dedup.bam > \
    ${SAMPLE}_mapping_depth.txt

# 3. Variant calling
# Create a bam index file
printf '\n[Step 9] Indexing bam file...\n'
bamtools index -in ${SAMPLE}_aln_sorted_dedup.bam

# Check BI/BD tags
# Lofreq needs special tags in BAM
# BI = insertion quality
# BD = deletion quality
# these represent confidence of indel observations.
# If tag exists:
  ## BI:Z:<<<<<<<<
  ## BD:Z:<<<<<<<<
printf '\n[Step 10] Check input BAM before running lofreq...\n'
if samtools view "${SAMPLE}_aln_sorted_dedup.bam" |
    grep -qm 1 -E '(^|[[:space:]])B[ID]:Z:'; then
    printf 'BI/BD tags found.\n'
else
    printf 'No BI/BD tags found; indel qualities may need to be added.\n'
fi

# Realign reads after removing duplicate reads
# this corrects mapping errors
printf '\n[Step 11] Realign reads in lofreq...\n'
lofreq viterbi \
  --verbose \
  -f ${REF}  \
  ${SAMPLE}_aln_sorted_dedup.bam | samtools sort > ${SAMPLE}_lofreq_realign.bam 
samtools index ${SAMPLE}_lofreq_realign.bam

# Run this if BI/BD tags are missing
# usually for bacteria/virus PE reads
# this inserts indel qualities 
# (=! indel 'alignment' qualities) just like GATK BQSR
printf '\n[Step 12] Insert indel qualities in lofreq...\n'
nice lofreq indelqual \
    --verbose \
    --dindel \
    -f ${REF} \
    -o ${SAMPLE}_lofreq.indelqual.bam \
    ${SAMPLE}_lofreq_realign.bam
samtools index ${SAMPLE}_lofreq.indelqual.bam

# Add base alignment quality (BAQ), incl indels
printf '\n[Step 13] Insert alignment quality in lofreq...\n'
lofreq alnqual -b \
    ${SAMPLE}_lofreq.indelqual.bam \
    ${REF} > ${SAMPLE}_lofreq_realign_alnq.bam
samtools index ${SAMPLE}_lofreq_realign_alnq.bam

# Variant calling
# this will add BAQ but not indel align qual
# thus recommend running lofreq alnqual beforehand.
printf '\n[Step 14] Calling variants in lofreq...\n'
nice lofreq call \
    --call-indels \
    --verbose \
    -f ${REF} \
    -o ${SAMPLE}_lofreq.vcf \
    ${SAMPLE}_lofreq_realign_alnq.bam
    # --no-default-filter
    # -a: set p-value threshold [default: 0.01 = Q20]


echo "LoFreq variant calling completed for ${SAMPLE}"
