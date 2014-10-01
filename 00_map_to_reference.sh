#!/bin/bash

set -e

# Build reference index:
/usr/local/src/bowtie2-2.2.3/bowtie2-build reference/Pc_ref.fa reference/Pc_ref > reference/bowtie2-build.out 2>&1 &
wait

# Map each sample to reference:
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc1_0d.fastq.gz -S mapped_reads/Pc1_0d.sam > mapped_reads/Pc1_0d.out 2>&1 &
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc1_1d.fastq.gz -S mapped_reads/Pc1_1d.sam > mapped_reads/Pc1_1d.out 2>&1 &
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc1_3d.fastq.gz -S mapped_reads/Pc1_3d.sam > mapped_reads/Pc1_3d.out 2>&1 &
wait
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc1_5d.fastq.gz -S mapped_reads/Pc1_5d.sam > mapped_reads/Pc1_5d.out 2>&1 &
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc1_9d.fastq.gz -S mapped_reads/Pc1_9d.sam > mapped_reads/Pc1_9d.out 2>&1 &
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc2_0d.fastq.gz -S mapped_reads/Pc2_0d.sam > mapped_reads/Pc2_0d.out 2>&1 &
wait
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc2_1d.fastq.gz -S mapped_reads/Pc2_1d.sam > mapped_reads/Pc2_1d.out 2>&1 &
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc2_3d.fastq.gz -S mapped_reads/Pc2_3d.sam > mapped_reads/Pc2_3d.out 2>&1 &
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc2_5d.fastq.gz -S mapped_reads/Pc2_5d.sam > mapped_reads/Pc2_5d.out 2>&1 &
wait
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc2_7d.fastq.gz -S mapped_reads/Pc2_7d.sam > mapped_reads/Pc2_7d.out 2>&1 &
/usr/local/src/bowtie2-2.2.3/bowtie2 --phred33 --very-sensitive-local -a -p 8 -x reference/Pc_ref -U data/Pc2_9d.fastq.gz -S mapped_reads/Pc2_9d.sam > mapped_reads/Pc2_9d.out 2>&1 &
wait
