#!/bin/bash

set -e

# Calculate read count per transcript by using Bowtie mapping.
python bowtie_map_to_counts.py mapped_reads/Pc1_0d.sam > read_counts/Pc1_0d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc1_1d.sam > read_counts/Pc1_1d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc1_3d.sam > read_counts/Pc1_3d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc1_5d.sam > read_counts/Pc1_5d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc1_9d.sam > read_counts/Pc1_9d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc2_0d.sam > read_counts/Pc2_0d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc2_1d.sam > read_counts/Pc2_1d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc2_3d.sam > read_counts/Pc2_3d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc2_5d.sam > read_counts/Pc2_5d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc2_7d.sam > read_counts/Pc2_7d.counts &
python bowtie_map_to_counts.py mapped_reads/Pc2_9d.sam > read_counts/Pc2_9d.counts &
wait
