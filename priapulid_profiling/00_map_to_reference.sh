#!/bin/bash

# Build reference index:
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2-build reference/Pc_ref.fa reference/Pc_ref > reference/bowtie2-build.log 2>&1 &

# Map each sample to reference:
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc1_0d.txt.gz -S mapped_reads/Pc1_0d_bowtie.map > mapped_reads/Pc1_0d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc1_1d.txt.gz -S mapped_reads/Pc1_1d_bowtie.map > mapped_reads/Pc1_1d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc1_3d.txt.gz -S mapped_reads/Pc1_3d_bowtie.map > mapped_reads/Pc1_3d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc1_5d.txt.gz -S mapped_reads/Pc1_5d_bowtie.map > mapped_reads/Pc1_5d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc1_9d.txt.gz -S mapped_reads/Pc1_9d_bowtie.map > mapped_reads/Pc1_9d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc2_0d.txt.gz -S mapped_reads/Pc2_0d_bowtie.map > mapped_reads/Pc2_0d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc2_1d.txt.gz -S mapped_reads/Pc2_1d_bowtie.map > mapped_reads/Pc2_1d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc2_3d.txt.gz -S mapped_reads/Pc2_3d_bowtie.map > mapped_reads/Pc2_3d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc2_5d.txt.gz -S mapped_reads/Pc2_5d_bowtie.map > mapped_reads/Pc2_5d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc2_7d.txt.gz -S mapped_reads/Pc2_7d_bowtie.map > mapped_reads/Pc2_7d.log 2>&1 &
/usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc2_9d.txt.gz -S mapped_reads/Pc2_9d_bowtie.map > mapped_reads/Pc2_9d.log 2>&1 &

