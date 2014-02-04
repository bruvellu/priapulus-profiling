Analyses
========

Based on [Helm et al. 2013](http://www.biomedcentral.com/1471-2164/14/266)

---

All 217.8 million reads from experiment (11 samples combined) and 464.4 million
reads from reference reads passed Illumina chastity filter. Reference assembled
with Agalma 0.3.5.

Mapping reads
-------------

Build reference index:

    /usr/local/src/bowtie2-2.0.0-beta7/bowtie2-build reference/Pc_ref.fa reference/Pc_ref > reference/bowtie2-build.log 2>&1 &

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

Calculate read counts
---------------------

Script to calculate read counts for each reference transcript. Slightly
modified from:
http://www.biomedcentral.com/content/supplementary/1471-2164-14-266-s8.py

```python
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
```

Ran with these commands:

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

Importing data to R
-------------------

Following steps from: http://www.biomedcentral.com/content/supplementary/1471-2164-14-266-s9.r

Create a table for each sample using read counts per transcript as input.

    Pc1_0d <- read.table("read_count_per_transcript/Pc1_0d.counts", header=TRUE)
    Pc1_1d <- read.table("read_count_per_transcript/Pc1_1d.counts", header=TRUE)
    Pc1_3d <- read.table("read_count_per_transcript/Pc1_3d.counts", header=TRUE)
    Pc1_5d <- read.table("read_count_per_transcript/Pc1_5d.counts", header=TRUE)
    Pc1_9d <- read.table("read_count_per_transcript/Pc1_9d.counts", header=TRUE)
    Pc2_0d <- read.table("read_count_per_transcript/Pc2_0d.counts", header=TRUE)
    Pc2_1d <- read.table("read_count_per_transcript/Pc2_1d.counts", header=TRUE)
    Pc2_3d <- read.table("read_count_per_transcript/Pc2_3d.counts", header=TRUE)
    Pc2_5d <- read.table("read_count_per_transcript/Pc2_5d.counts", header=TRUE)
    Pc2_7d <- read.table("read_count_per_transcript/Pc2_7d.counts", header=TRUE)
    Pc2_9d <- read.table("read_count_per_transcript/Pc2_9d.counts", header=TRUE)

Merge data into matrix by common row names (transcripts).

    merged_data <- merge(Pc1_0d, Pc2_0d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d")
    merged_data <- merge(merged_data, Pc1_1d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d")
    merged_data <- merge(merged_data, Pc2_1d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d")
    merged_data <- merge(merged_data, Pc1_3d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d")
    merged_data <- merge(merged_data, Pc2_3d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d", "count_Pc2_3d")
    merged_data <- merge(merged_data, Pc1_5d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d", "count_Pc2_3d", "count_Pc1_5d")
    merged_data <- merge(merged_data, Pc2_5d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d", "count_Pc2_3d", "count_Pc1_5d", "count_Pc2_5d")
    merged_data <- merge(merged_data, Pc2_7d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d", "count_Pc2_3d", "count_Pc1_5d", "count_Pc2_5d", "count_Pc2_7d")
    merged_data <- merge(merged_data, Pc1_9d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d", "count_Pc2_3d", "count_Pc1_5d", "count_Pc2_5d", "count_Pc2_7d", "count_Pc1_9d")
    merged_data <- merge(merged_data, Pc2_9d, by.x="reference", by.y="reference", all=TRUE)
    names(merged_data) <- c("transcript", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d", "count_Pc2_3d", "count_Pc1_5d", "count_Pc2_5d", "count_Pc2_7d", "count_Pc1_9d", "count_Pc2_9d")

Name rows with transcript name.

    rownames(merged_data) <- merged_data[,1]

Keep only numeric columns.

    n <- data.matrix(merged_data[,2:12])

Set NA values to 0:

    n[is.na(n)] <- 0

Prepare data with edgeR
-----------------------

Calculate normalization factor for each column. Read more about normalization methods [here](http://www.biomedcentral.com/1471-2105/11/94).

    # TODO Which normalization method is the default???
    nf <- calcNormFactors(n)

Get sum of each column.

    lib.size <- colSums(n)

Get effective size by multiplying size by normalization factor.

    lib.effective.size <- lib.size * nf

Get normalization multiplier to be applied into original counts.

    norm.multiplier <- 1000000 / lib.effective.size

Create normalized matrix by multiplying normalization factor to counts.

    q <- n * norm.multiplier

Rename columns.

    colnames(q) <- c("norm_count_Pc1_0d", "norm_count_Pc2_0d", "norm_count_Pc1_1d", "norm_count_Pc2_1d", "norm_count_Pc1_3d", "norm_count_Pc2_3d", "norm_count_Pc1_5d", "norm_count_Pc2_5d", "norm_count_Pc2_7d", "norm_count_Pc1_9d", "norm_count_Pc2_9d")

Bind counts and normalize data.

    merged_data_norm <- cbind(merged_data, q)

Remove NAs again.

    merged_data_norm[is.na(merged_data_norm)] <- 0

STEM analysis
-------------

Calculate the average between replicates.

    avg_norm_count_0d <- (merged_data_norm$norm_count_Pc1_0d + merged_data_norm$norm_count_Pc2_0d) / 2
    avg_norm_count_1d <- (merged_data_norm$norm_count_Pc1_1d + merged_data_norm$norm_count_Pc2_1d) / 2
    avg_norm_count_3d <- (merged_data_norm$norm_count_Pc1_3d + merged_data_norm$norm_count_Pc2_3d) / 2
    avg_norm_count_5d <- (merged_data_norm$norm_count_Pc1_5d + merged_data_norm$norm_count_Pc2_5d) / 2
    avg_norm_count_7d <- merged_data_norm$norm_count_Pc2_7d
    avg_norm_count_9d <- (merged_data_norm$norm_count_Pc1_9d + merged_data_norm$norm_count_Pc2_9d) / 2

Build data frame for average values.

    avg_merged_data_norm <- cbind(avg_norm_count_0d, avg_norm_count_1d, avg_norm_count_3d, avg_norm_count_5d, avg_norm_count_7d, avg_norm_count_9d)

Put names on rows.

    rownames(avg_merged_data_norm) <- merged_data_norm[,1]

Write file with average data for STEM input.

    write.table(avg_merged_data_norm, "average")

Edit the file `average` to run STEM.

    cp average avg_stem_input
    vim avg_stem_input

1. Manually add `'"transcript" '` as the first column name (without single quotes, note the white space).
2. Remove quotes with `:%s:"::g`.
3. Substitute white space for tabs `:%s:\s\+:\t:g`.

Ran STEM with the following command:

    java -mx1024M -jar ~/src/stem/stem.jar

Default settings used as shown below (see [complete output](stem/stem_output)):

    #Main Input:
    Data_File   /home/nelas/Biologia/Doutorado/Priapulus/rna/RNAseq_profiling/stem/avg_stem_input
    Gene_Annotation_Source  No annotations
    Gene_Annotation_File
    Cross_Reference_Source  No cross references
    Cross_Reference_File
    Gene_Location_Source    No Gene Locations
    Gene_Location_File
    Clustering_Method[STEM Clustering Method,K-means]   STEM Clustering Method
    Maximum_Number_of_Model_Profiles    50
    Maximum_Unit_Change_in_Model_Profiles_between_Time_Points   2
    Normalize_Data[Log normalize data,Normalize data,No normalization/add 0]    Normalize data
    Spot_IDs_included_in_the_data_file  false

From the initial 58133 transcripts, STEM filtered out 33023 while 25110 passed
the filter. Below are the profiles found ordered by transcript abundance and
significance:

![Priapulus caudatus STEM Profiles](stem/Pcau_stem.png)

Relevant profiles:

| profile                                        | description                               |
| :------:                                       | :----------                               |
| [profile 8](profile_8)                         | constant decrease from initial oocytes.   |
| [profile 39](profile_39)                       | constant increase from initial oocytes.   |
| [profile 31](profile_31)                       | peak at 3d coincident with gastrulation.  |
| [profile 25](profile_25)                       | peak at 3d but drastically downregulated. |
| [profile 22](profile_22)                       | peak at 7d coincident with introvertula.  |
| [profiles 17](profile_17) and [18](profile_18) | low during cleavage, then up.             |
| [profile 1](profile_1)                         | low expression since cleavage.            |

[profile 8]: stem/profile_8
[profile 39]: stem/profile_39
[profile 31]: stem/profile_31
[profile 25]: stem/profile_25
[profile 22]: stem/profile_22
[profiles 17]: stem/profile_17
[profile_18]: stem/profile_18
[profile 1]: stem/profile_1
