Analyses
========

Based on [Helm et al. 2013](http://www.biomedcentral.com/1471-2164/14/266)

---

All 217.8 million reads from experiment (11 samples combined) and 464.4 million
reads from reference reads passed Illumina chastity filter. Reference assembled
with Agalma 0.3.5.

Mapping reads
-------------

Build reference index using `bowtie2-build`:

    /usr/local/src/bowtie2-2.0.0-beta7/bowtie2-build reference/Pc_ref.fa reference/Pc_ref > reference/bowtie2-build.log 2>&1 &

Map each sample to reference. Example below:

    /usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc1_0d.txt.gz -S mapped_reads/Pc1_0d_bowtie.map > mapped_reads/Pc1_0d.log 2>&1 &
    ...

**Source:** [00_map_to_reference.sh](00_map_to_ref.sh)

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

**Source:** [bowtie_map_to_counts.py](bowtie_map_to_counts.py)

Ran command for each sample:

    ./bowtie_map_to_counts.py mapped_reads/Pc1_0d_bowtie.map > read_count_per_transcript/Pc1_0d.counts &
    ...

**Source:** [01_calculate_read_count.sh](01_calculate_read_count.sh)

Importing data to R
-------------------

Following steps from: http://www.biomedcentral.com/content/supplementary/1471-2164-14-266-s9.r

```s
# Load read counts.
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

# Merge read counts intersecting by transcript.
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
names(merged_data) <- c("reference", "count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d", "count_Pc2_3d", "count_Pc1_5d", "count_Pc2_5d", "count_Pc2_7d", "count_Pc1_9d", "count_Pc2_9d")

# Name rows with transcript name.
rownames(merged_data) <- merged_data[,1]
```

**Source:** [02_import_to_r.sh](02_import_to_r.sh) and [load_read_counts.r](load_read_counts.r)

Normalize data with edgeR
-------------------------

```s
# Keep only numeric columns.
num_data <- data.matrix(merged_data[,2:12])

# Set NA values to 0:
num_data[is.na(num_data)] <- 0

# Calculate normalization factor for each column.
norm_factors <- calcNormFactors(num_data)

# Get sum of each column.
lib.size <- colSums(num_data)

# Get effective size by multiplying size by normalization factor.
lib.effective.size <- lib.size * norm_factors

# Get normalization multiplier to be applied into original counts.
norm.multiplier <- 1000000 / lib.effective.size

# Create normalized matrix by multiplying normalization factor to counts.
norm_data <- num_data * norm.multiplier

# Rename columns.
colnames(norm_data) <- c("norm_count_Pc1_0d", "norm_count_Pc2_0d", "norm_count_Pc1_1d", "norm_count_Pc2_1d", "norm_count_Pc1_3d", "norm_count_Pc2_3d", "norm_count_Pc1_5d", "norm_count_Pc2_5d", "norm_count_Pc2_7d", "norm_count_Pc1_9d", "norm_count_Pc2_9d")

# Bind counts and normalize data.
merged_data <- cbind(merged_data, norm_data)

# Remove NAs again.
merged_data[is.na(merged_data)] <- 0
```

**Source:** [03_normalize.sh](03_normalize.sh) and [edger_normalize.r](edger_normalize.r)

STEM analysis
-------------

```s
# Calculate the average between replicates.
avg_norm_count_0d <- (merged_data$norm_count_Pc1_0d + merged_data$norm_count_Pc2_0d) / 2
avg_norm_count_1d <- (merged_data$norm_count_Pc1_1d + merged_data$norm_count_Pc2_1d) / 2
avg_norm_count_3d <- (merged_data$norm_count_Pc1_3d + merged_data$norm_count_Pc2_3d) / 2
avg_norm_count_5d <- (merged_data$norm_count_Pc1_5d + merged_data$norm_count_Pc2_5d) / 2
avg_norm_count_7d <- merged_data$norm_count_Pc2_7d
avg_norm_count_9d <- (merged_data$norm_count_Pc1_9d + merged_data$norm_count_Pc2_9d) / 2

# Build data frame for average values.
avg_data_norm <- cbind(avg_norm_count_0d, avg_norm_count_1d, avg_norm_count_3d, avg_norm_count_5d, avg_norm_count_7d, avg_norm_count_9d)

# Put names on rows.
rownames(avg_data_norm) <- merged_data[,1]

# Bind average counts with matrix.
merged_data <- cbind(merged_data, avg_data_norm)

# Write file with average data for STEM input.
write.table(avg_data_norm, "average")
```

**Source:** [04_stem_profiles.sh](04_stem_profiles.sh) and [stem_prepare.r](stem_prepare.r)

Edit the file `average` to run STEM.

    cp average avg_stem_input
    vim avg_stem_input

1. Manually add `'"transcript" '` as the first column name (without single quotes, note the white space).
2. Remove quotes with `:%s:"::g`.
3. Substitute white space for tabs `:%s:\s\+:\t:g`.

Run STEM with the following command:

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
| [profile 8][profile_8]                         | constant decrease from initial oocytes    |
| [profile 39][profile_39]                       | constant increase from initial oocytes    |
| [profile 31][profile_31]                       | peak at 3d coincident with gastrulation   |
| [profile 25][profile_25]                       | peak at 3d but drastically downregulated  |
| [profile 22][profile_22]                       | peak at 7d coincident with introvertula   |
| [profile 17][profile_17] and [18][profile_18]  | low during cleavage, then up              |
| [profile 1][profile_1]                         | low expression since cleavage             |

[profile_8]: stem/profile_8
[profile_39]: stem/profile_39
[profile_31]: stem/profile_31
[profile_25]: stem/profile_25
[profile_22]: stem/profile_22
[profile_17]: stem/profile_17
[profile_18]: stem/profile_18
[profile_1]: stem/profile_1

Output file was created manually by saving the Main Table for Genes Passing
Filter as `genes_passing_filter`. This file needs to be edited to be imported
back into R.

    cp stem_output stem_profiles_to_r

1. Remove tab from beginning of the line with `:%s:^\t::g`
2. Remove initial header with `:%s:^Selected\t::g`
3. Remove SPOT column with `%s:\tID_\d\+\t\(\d\+\):\t\1:g`
4. Remove SPOT header with `:%s:SPOT\t::g`
5. Change count to delta in header: `:%s:avg_norm_count:avg_norm_delta:g`
6. Remove first line: `Table of Genes Passing Filter`
7. Lowercase everything with: `ggVGu`
8. Uppercase locus and transcript on deflines with `:%s:locus_:Locus_:g` and `:%s:transcript_:Transcript_:g`

```s
# Put STEM data into matrix with R
stem <- read.table("stem_plot_data", header = TRUE)

Zlower <- tolower(Z$transcript)
Zlower <- as.data.frame(Zlower)
rownames(Z) <- Zlower[,1]
colnames(Zs)[1] <- "lowercase_transcript"
Zs <- merge(Z, stem, by.x= "row.names", by.y = "lowercase_transcript", all=TRUE)

#lowercase_transcript was a temporary column to merge STEM data.

# Saving full database in R:

write.table(Zs, "matrix.R")
```

Replicate plots
---------------

Plotting the normalized average count between replicates, except for 7d sample.

![Scatter plots between replicates](plots/rep_plots.png)

**Source:** [05_plot_replicates.sh](05_plot_replicates.sh) and [build_scatter_plots.r](build_scatter_plots.r)

Differential expression
-----------------------

Gene ontology
-------------

1. Parse accession numbers from annotated assembly.
2. Use Entrez to fetch list of GO ids associated with the protein. (should I fetch InterPro?)
3. For each GO id fetch hierarchy, terms, and evidence codes.
4. Associate terms with contigs.
5. Done.
