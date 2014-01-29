Analyses
========

Based on Helm et al. 2013

---

All 217.8 million reads from experiment (11 samples combined) and 464.4 million
reads from reference reads passed Illumina chastity filter. Reference assembled
with Agalma 0.3.5.

Mapping reads
-------------

Build reference index:

    /usr/local/src/bowtie2-2.0.0-beta7/bowtie2-build reference/Pc_ref.fa reference/Pc_ref

Map each sample to reference:

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

<!--Attach report?-->

Calculate read counts
---------------------

Script to calculate read counts for each reference transcript. Slightly
modified from:
http://www.biomedcentral.com/content/supplementary/1471-2164-14-266-s8.py

    from collections import defaultdict

    map_file = open('mapped_reads/Pc1_1d_bowtie.map')
    map = defaultdict(set)

    for line in map_file:
        if not line.startswith('@'):
            line = line.strip()
            fields = line.split('\t')
            read = fields[0]
            ref = fields[2]
            map[read].add(ref)

    counts = defaultdict(int)

    for key, hits in map.iteritems():
        if len(hits) > 1:
            continue
        ref = list(hits)[0]
        counts[ref] = counts[ref] + 1

    print "reference\tcount"
    for ref, count in counts.iteritems():
        print ref + "\t" + str(count)

Ran with these commands:

    ./bowtie_map_to_counts.py ../mapped_reads/Pc1_0d_bowtie.map > Pc1_0d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc1_1d_bowtie.map > Pc1_1d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc1_3d_bowtie.map > Pc1_3d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc1_5d_bowtie.map > Pc1_5d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc1_9d_bowtie.map > Pc1_9d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc2_0d_bowtie.map > Pc2_0d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc2_1d_bowtie.map > Pc2_1d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc2_3d_bowtie.map > Pc2_3d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc2_5d_bowtie.map > Pc2_5d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc2_7d_bowtie.map > Pc2_7d.counts &
    ./bowtie_map_to_counts.py ../mapped_reads/Pc2_9d_bowtie.map > Pc2_9d.counts &

