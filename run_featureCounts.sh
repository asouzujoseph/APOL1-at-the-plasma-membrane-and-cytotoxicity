#!/bin/bash

# Load mamba and env
source ~/.local/share/mamba/etc/profile.d/mamba.sh
mamba activate rnaseq

# Directories
ALIGN=/mnt/c/Users/User/Documents/Damola/aligned
GTF=/mnt/c/Users/User/Documents/Damola/references/Homo_sapiens/UCSC/hg38/Annotation/GENCODE/gencode.v49.annotation.gtf
OUT=/mnt/c/Users/User/Documents/Damola/counts

mkdir -p "$OUT"

echo "Running featureCounts on all BAM files in: $ALIGN"

featureCounts \
  -T 4 \
  -p \
  -B \
  -C \
  -a "$GTF" \
  -o "$OUT/gene_counts.txt" \
  "$ALIGN"/*.bam

echo "featureCounts complete. Output: $OUT/gene_counts.txt"

mamba deactivate
