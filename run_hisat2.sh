#!/bin/bash

source ~/.local/share/mamba/etc/profile.d/mamba.sh
mamba activate rnaseq

TRIM=/mnt/c/Users/User/Documents/Damola/trimmed
ALIGN=/mnt/c/Users/User/Documents/Damola/aligned
INDEX=/mnt/c/Users/User/Documents/Damola/references/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/hisat2_index/genome

mkdir -p "$ALIGN"

echo "Running HISAT2 alignment using 1 job (2 threads each)..."

ls $TRIM/*_1_trimmed.fastq.gz | parallel -j 1 '
    R1={};
    R2=${R1/_1_trimmed.fastq.gz/_2_trimmed.fastq.gz};
    base=$(basename $R1 _1_trimmed.fastq.gz);

    hisat2 -p 2 -x '"$INDEX"' \
        -1 $R1 -2 $R2 \
        | samtools sort -@ 2 -o '"$ALIGN"'/${base}.bam
'

echo "HISAT2 alignment complete. BAM files saved in: $ALIGN"

mamba deactivate
