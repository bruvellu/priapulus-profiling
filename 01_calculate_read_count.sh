# Calculate read count per transcript by using Bowtie mapping.
./bowtie_map_to_counts.py mapped_reads/Pc1_0d_bowtie.map > read_count_per_transcript/Pc1_0d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc1_1d_bowtie.map > read_count_per_transcript/Pc1_1d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc1_3d_bowtie.map > read_count_per_transcript/Pc1_3d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc1_5d_bowtie.map > read_count_per_transcript/Pc1_5d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc1_9d_bowtie.map > read_count_per_transcript/Pc1_9d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc2_0d_bowtie.map > read_count_per_transcript/Pc2_0d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc2_1d_bowtie.map > read_count_per_transcript/Pc2_1d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc2_3d_bowtie.map > read_count_per_transcript/Pc2_3d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc2_5d_bowtie.map > read_count_per_transcript/Pc2_5d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc2_7d_bowtie.map > read_count_per_transcript/Pc2_7d.counts &
./bowtie_map_to_counts.py mapped_reads/Pc2_9d_bowtie.map > read_count_per_transcript/Pc2_9d.counts &

