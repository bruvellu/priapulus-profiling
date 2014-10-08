#!/bin/bash

# $1 = SAM file
# $2 = BAM output without extension

SAM="$1"
BAM="$2.bam"
BAM_SORTED="$2_sorted"
BAI="$2_sorted.bai"

# Generate BAM from SAM.
samtools view -b -S $SAM > $BAM

# Sort BAM file.
samtools sort $BAM $BAM_SORTED

# Build index for BAM.
samtools index "$BAM_SORTED.bam" $BAI
