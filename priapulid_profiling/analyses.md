Analyses
========

Based on [Helm et al. 2013](http://www.biomedcentral.com/1471-2164/14/266)

---

All 217.8 million reads from experiment (11 samples combined) and 464.4 million
reads from reference reads passed Illumina chastity filter. Reference assembled
with Agalma 0.3.5.


Summary
-------

1. [00_map_to_reference.sh](00_map_to_ref.sh) Build reference assembly index
and map reads to it with Bowtie2.
2. [01_calculate_read_count.sh](01_calculate_read_count.sh) Calculate read
counts for each transcript of every sample.

Mapping reads
-------------

Build reference index using `bowtie2-build`:

    /usr/local/src/bowtie2-2.0.0-beta7/bowtie2-build reference/Pc_ref.fa reference/Pc_ref > reference/bowtie2-build.log 2>&1 &

Map each sample to reference. Example below:

    /usr/local/src/bowtie2-2.0.0-beta7/bowtie2 --phred33 --very-sensitive-local -a -p 10 -x reference/Pc_ref -U data/Pc1_0d.txt.gz -S mapped_reads/Pc1_0d_bowtie.map > mapped_reads/Pc1_0d.log 2>&1 &
    ...

Calculate read counts
---------------------

Calculate the number of reads mapping to each reference transcript. I used a
slightly modified version of [Helm et al. 2013
bowtie_map_to_counts][Helm2013_bowtie]

**Source:** [bowtie_map_to_counts.py](bowtie_map_to_counts.py)

Run command for each sample:

    ./bowtie_map_to_counts.py mapped_reads/Pc1_0d_bowtie.map > read_count_per_transcript/Pc1_0d.counts &
    ...

[Helm2013_bowtie]: http://www.biomedcentral.com/content/supplementary/1471-2164-14-266-s8.py

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

# TODO Remember to remove counts for unmapped transcripts (*).
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
stem <- read.table("stem_profiles_to_r", header=TRUE)
merged_data <- merge(merged_data, stem, by.x="reference", by.y="transcript", all=TRUE)

# Save matrix in R:
write.table(merged_data, "matrix.r")
```

Replicate plots
---------------

Plotting the normalized average count between replicates, except for 7d sample.

![Scatter plots between replicates](plots/rep_plots.png)

**Source:** [05_plot_replicates.sh](05_plot_replicates.sh) and [build_scatter_plots.r](build_scatter_plots.r)

Differential expression
-----------------------

This scripts runs step-by-step the test for differential gene expression as follows:

1. Subset read counts by interval and define available time points.
2. Create DGEList object for running the test.
3. Calculate the normalization factor and dispersion.
4. Run test and output a table with counts of significant genes.

```s
# Load edgeR.
library(edgeR)

# TODO write single function to do the exactTest and pass on arguments for each
# experimental time points interval.

## 0d x 1d

# Subset counts from matrix for 0d-1d test.
data_0d1d <- merged_data[c("count_Pc1_0d", "count_Pc2_0d", "count_Pc1_1d", "count_Pc2_1d")]
rownames(data_0d1d) <- merged_data$reference

# Set time points.
time_points <- factor(c("0d", "0d", "1d", "1d"))
tp_0d1d <- data.frame(Sample=colnames(data_0d1d), time_points)

# Object for running the test.
dge_0d1d <- DGEList(counts=data_0d1d, group=tp_0d1d$time_points)
colnames(dge_0d1d) <- tp_0d1d$Sample

# Filter counts?
##filter read with no/low counts
#keep <- rowSums(cpm(y)) >= 2
#y <- y[keep,]
#
##update library sizes
#y$samples$lib.size <- colSums(y$counts)
#dim(y)

# Calculate normalization factor.
dge_0d1d <- calcNormFactors(dge_0d1d)

# Estimate dispersion.
dge_0d1d <- estimateCommonDisp(dge_0d1d, verbose=TRUE)
dge_0d1d <- estimateTagwiseDisp(dge_0d1d)

# Test for differential gene expression.
test_0d1d <- exactTest(dge_0d1d, pair=c("0d", "1d"))

# Build table of results.
table_0d1d <- test_0d1d$table
table_0d1d <- cbind(rownames(table_0d1d), table_0d1d)
names(table_0d1d)[names(table_0d1d) == "rownames(table_0d1d)"] = "transcript"

# Adjust p-values with bonferroni.
table_0d1d <- cbind(table_0d1d, p.adjust(table_0d1d$PValue, method="bonferroni"))
names(table_0d1d)[names(table_0d1d)=="p.adjust(table_0d1d$PValue, method = \"bonferroni\")"] = "PAdjust_0d_1d"
names(table_0d1d)[names(table_0d1d)=="logFC"] = "logFC_0d_1d"
names(table_0d1d)[names(table_0d1d)=="logCPM"] = "logCPM_0d_1d"
names(table_0d1d)[names(table_0d1d)=="PValue"] = "PValue_0d_1d"

## 1d x 3d

# Subset counts from matrix for 1d-3d test.
data_1d3d <- merged_data[c("count_Pc1_1d", "count_Pc2_1d", "count_Pc1_3d", "count_Pc2_3d")]
rownames(data_1d3d) <- merged_data$reference

# Set time points.
time_points <- factor(c("1d", "1d", "3d", "3d"))
tp_1d3d <- data.frame(Sample=colnames(data_1d3d), time_points)

# Object for running the test.
dge_1d3d <- DGEList(counts=data_1d3d, group=tp_1d3d$time_points)
colnames(dge_1d3d) <- tp_1d3d$Sample

# Filter counts?
##filter read with no/low counts
#keep <- rowSums(cpm(y)) >= 2
#y <- y[keep,]
#
##update library sizes
#y$samples$lib.size <- colSums(y$counts)
#dim(y)

# Calculate normalization factor.
dge_1d3d <- calcNormFactors(dge_1d3d)

# Estimate dispersion.
dge_1d3d <- estimateCommonDisp(dge_1d3d, verbose=TRUE)
dge_1d3d <- estimateTagwiseDisp(dge_1d3d)

# Test for differential gene expression.
test_1d3d <- exactTest(dge_1d3d, pair=c("1d", "3d"))

# Build table of results.
table_1d3d <- test_1d3d$table
table_1d3d <- cbind(rownames(table_1d3d), table_1d3d)
names(table_1d3d)[names(table_1d3d) == "rownames(table_1d3d)"] = "transcript"

# Adjust p-values with bonferroni.
table_1d3d <- cbind(table_1d3d, p.adjust(table_1d3d$PValue, method="bonferroni"))
names(table_1d3d)[names(table_1d3d)=="p.adjust(table_1d3d$PValue, method = \"bonferroni\")"] = "PAdjust_1d_3d"
names(table_1d3d)[names(table_1d3d)=="logFC"] = "logFC_1d_3d"
names(table_1d3d)[names(table_1d3d)=="logCPM"] = "logCPM_1d_3d"
names(table_1d3d)[names(table_1d3d)=="PValue"] = "PValue_1d_3d"

## 3d x 5d

# Subset counts from matrix for 3d-5d test.
data_3d5d <- merged_data[c("count_Pc1_3d", "count_Pc2_3d", "count_Pc1_5d", "count_Pc2_5d")]
rownames(data_3d5d) <- merged_data$reference

# Set time points.
time_points <- factor(c("3d", "3d", "5d", "5d"))
tp_3d5d <- data.frame(Sample=colnames(data_3d5d), time_points)

# Object for running the test.
dge_3d5d <- DGEList(counts=data_3d5d, group=tp_3d5d$time_points)
colnames(dge_3d5d) <- tp_3d5d$Sample

# Filter counts?
##filter read with no/low counts
#keep <- rowSums(cpm(y)) >= 2
#y <- y[keep,]
#
##update library sizes
#y$samples$lib.size <- colSums(y$counts)
#dim(y)

# Calculate normalization factor.
dge_3d5d <- calcNormFactors(dge_3d5d)

# Estimate dispersion.
dge_3d5d <- estimateCommonDisp(dge_3d5d, verbose=TRUE)
dge_3d5d <- estimateTagwiseDisp(dge_3d5d)

# Test for differential gene expression.
test_3d5d <- exactTest(dge_3d5d, pair=c("3d", "5d"))

# Build table of results.
table_3d5d <- test_3d5d$table
table_3d5d <- cbind(rownames(table_3d5d), table_3d5d)
names(table_3d5d)[names(table_3d5d) == "rownames(table_3d5d)"] = "transcript"

# Adjust p-values with bonferroni.
table_3d5d <- cbind(table_3d5d, p.adjust(table_3d5d$PValue, method="bonferroni"))
names(table_3d5d)[names(table_3d5d)=="p.adjust(table_3d5d$PValue, method = \"bonferroni\")"] = "PAdjust_3d_5d"
names(table_3d5d)[names(table_3d5d)=="logFC"] = "logFC_3d_5d"
names(table_3d5d)[names(table_3d5d)=="logCPM"] = "logCPM_3d_5d"
names(table_3d5d)[names(table_3d5d)=="PValue"] = "PValue_3d_5d"

## 5d x 7d

# Subset counts from matrix for 5d-7d test.
data_5d7d <- merged_data[c("count_Pc1_5d", "count_Pc2_5d", "count_Pc2_7d")]
rownames(data_5d7d) <- merged_data$reference

# Set time points.
time_points <- factor(c("5d", "5d", "7d"))
tp_5d7d <- data.frame(Sample=colnames(data_5d7d), time_points)

# Object for running the test.
dge_5d7d <- DGEList(counts=data_5d7d, group=tp_5d7d$time_points)
colnames(dge_5d7d) <- tp_5d7d$Sample

# Filter counts?
##filter read with no/low counts
#keep <- rowSums(cpm(y)) >= 2
#y <- y[keep,]
#
##update library sizes
#y$samples$lib.size <- colSums(y$counts)
#dim(y)

# Calculate normalization factor.
dge_5d7d <- calcNormFactors(dge_5d7d)

# Estimate dispersion.
dge_5d7d <- estimateCommonDisp(dge_5d7d, verbose=TRUE)
dge_5d7d <- estimateTagwiseDisp(dge_5d7d)

# Test for differential gene expression.
test_5d7d <- exactTest(dge_5d7d, pair=c("5d", "7d"))

# Build table of results.
table_5d7d <- test_5d7d$table
table_5d7d <- cbind(rownames(table_5d7d), table_5d7d)
names(table_5d7d)[names(table_5d7d) == "rownames(table_5d7d)"] = "transcript"

# Adjust p-values with bonferroni.
table_5d7d <- cbind(table_5d7d, p.adjust(table_5d7d$PValue, method="bonferroni"))
names(table_5d7d)[names(table_5d7d)=="p.adjust(table_5d7d$PValue, method = \"bonferroni\")"] = "PAdjust_5d_7d"
names(table_5d7d)[names(table_5d7d)=="logFC"] = "logFC_5d_7d"
names(table_5d7d)[names(table_5d7d)=="logCPM"] = "logCPM_5d_7d"
names(table_5d7d)[names(table_5d7d)=="PValue"] = "PValue_5d_7d"

## 7d x 9d

# Subset counts from matrix for 7d-9d test.
data_7d9d <- merged_data[c("count_Pc2_7d", "count_Pc1_9d", "count_Pc2_9d")]
rownames(data_7d9d) <- merged_data$reference

# Set time points.
time_points <- factor(c("7d", "9d", "9d"))
tp_7d9d <- data.frame(Sample=colnames(data_7d9d), time_points)

# Object for running the test.
dge_7d9d <- DGEList(counts=data_7d9d, group=tp_7d9d$time_points)
colnames(dge_7d9d) <- tp_7d9d$Sample

# Filter counts?
##filter read with no/low counts
#keep <- rowSums(cpm(y)) >= 2
#y <- y[keep,]
#
##update library sizes
#y$samples$lib.size <- colSums(y$counts)
#dim(y)

# Calculate normalization factor.
dge_7d9d <- calcNormFactors(dge_7d9d)

# Estimate dispersion.
dge_7d9d <- estimateCommonDisp(dge_7d9d, verbose=TRUE)
dge_7d9d <- estimateTagwiseDisp(dge_7d9d)

# Test for differential gene expression.
test_7d9d <- exactTest(dge_7d9d, pair=c("7d", "9d"))

# Build table of results.
table_7d9d <- test_7d9d$table
table_7d9d <- cbind(rownames(table_7d9d), table_7d9d)
names(table_7d9d)[names(table_7d9d) == "rownames(table_7d9d)"] = "transcript"

# Adjust p-values with bonferroni.
table_7d9d <- cbind(table_7d9d, p.adjust(table_7d9d$PValue, method="bonferroni"))
names(table_7d9d)[names(table_7d9d)=="p.adjust(table_7d9d$PValue, method = \"bonferroni\")"] = "PAdjust_7d_9d"
names(table_7d9d)[names(table_7d9d)=="logFC"] = "logFC_7d_9d"
names(table_7d9d)[names(table_7d9d)=="logCPM"] = "logCPM_7d_9d"
names(table_7d9d)[names(table_7d9d)=="PValue"] = "PValue_7d_9d"

# Merge results from DGE.
merged_dge <- merge(table_0d1d, table_1d3d, by.x="transcript", by.y="transcript", all=TRUE)
merged_dge <- merge(merged_dge, table_3d5d, by.x="transcript", by.y="transcript", all=TRUE)
merged_dge <- merge(merged_dge, table_5d7d, by.x="transcript", by.y="transcript", all=TRUE)
merged_dge <- merge(merged_dge, table_7d9d, by.x="transcript", by.y="transcript", all=TRUE)

# Exclude previous edgeR values (if it is a re-run).
merged_data$PValue_0d_1d <- NULL
merged_data$PAdjust_0d_1d <- NULL
merged_data$logFC_0d_1d <- NULL
merged_data$logCPM_0d_1d <- NULL
merged_data$PValue_1d_3d <- NULL
merged_data$PAdjust_1d_3d <- NULL
merged_data$logFC_1d_3d <- NULL
merged_data$logCPM_1d_3d <- NULL
merged_data$PValue_3d_5d <- NULL
merged_data$PAdjust_3d_5d <- NULL
merged_data$logFC_3d_5d <- NULL
merged_data$logCPM_3d_5d <- NULL
merged_data$PValue_5d_7d <- NULL
merged_data$PAdjust_5d_7d <- NULL
merged_data$logFC_5d_7d <- NULL
merged_data$logCPM_5d_7d <- NULL
merged_data$PValue_7d_9d <- NULL
merged_data$PAdjust_7d_9d <- NULL
merged_data$logFC_7d_9d <- NULL
merged_data$logCPM_7d_9d <- NULL

# Merge differential expression into matrix.
merged_data <- merge(merged_data, merged_dge, by.x="reference", by.y="transcript", all=TRUE)
# TODO Is there a need to sort?

# Select differentially expressed genes for each interval.

# 0d - 1d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
i1 <- merged_data[(merged_data$logFC_0d_1d > 0) & (merged_data$PAdjust_0d_1d < 0.05) & (!is.na(merged_data$PAdjust_0d_1d)),]
d1 <- merged_data[(merged_data$logFC_0d_1d < 0) & (merged_data$PAdjust_0d_1d < 0.05) & (!is.na(merged_data$PAdjust_0d_1d)),]

# 1d - 3d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
i2 <- merged_data[(merged_data$logFC_1d_3d > 0) & (merged_data$PAdjust_1d_3d < 0.05) & (!is.na(merged_data$PAdjust_1d_3d)),]
d2 <- merged_data[(merged_data$logFC_1d_3d < 0) & (merged_data$PAdjust_1d_3d < 0.05) & (!is.na(merged_data$PAdjust_1d_3d)),]

# 3d - 5d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
i3 <- merged_data[(merged_data$logFC_3d_5d > 0) & (merged_data$PAdjust_3d_5d < 0.05) & (!is.na(merged_data$PAdjust_3d_5d)),]
d3 <- merged_data[(merged_data$logFC_3d_5d < 0) & (merged_data$PAdjust_3d_5d < 0.05) & (!is.na(merged_data$PAdjust_3d_5d)),]

# 5d - 7d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
i4 <- merged_data[(merged_data$logFC_5d_7d > 0) & (merged_data$PAdjust_5d_7d < 0.05) & (!is.na(merged_data$PAdjust_5d_7d)),]
d4 <- merged_data[(merged_data$logFC_5d_7d < 0) & (merged_data$PAdjust_5d_7d < 0.05) & (!is.na(merged_data$PAdjust_5d_7d)),]

# 7d - 9d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
i5 <- merged_data[(merged_data$logFC_7d_9d > 0) & (merged_data$PAdjust_7d_9d < 0.05) & (!is.na(merged_data$PAdjust_7d_9d)),]
d5 <- merged_data[(merged_data$logFC_7d_9d < 0) & (merged_data$PAdjust_7d_9d < 0.05) & (!is.na(merged_data$PAdjust_7d_9d)),]

# Count increasing and decreasing genes per interval.
I1 <- c(nrow(i1),-(nrow(d1)))
I2 <- c(nrow(i2),-(nrow(d2)))
I3 <- c(nrow(i3),-(nrow(d3)))
I4 <- c(nrow(i4),-(nrow(d4)))
I5 <- c(nrow(i5),-(nrow(d5)))

# Save counts into a data frame.
de_counts <- data.frame(I1, I2, I3, I4, I5)
```

Dispersion values for each interval:

| interval | dispersion | bcv
| :------: | :--------: | :-:
| 0d_1d    | 0.43092    | 0.6564
| 1d_3d    | 0.32958    | 0.5741
| 3d_5d    | 0.26231    | 0.5122
| 5d_7d    | 0.21773    | 0.4666
| 7d_9d    | 0.15213    | 0.39

**Source:** [test_dge.r](test_dge.r)

Plots for differentially expressed genes
----------------------------------------

![de_genes](de_genes.png)

```s
# Subset each interval keeping only valid p-values.
m1 <- subset(merged_data, !is.na(merged_data$PAdjust_0d_1d))
rownames(m1)=m1[,1]
m2 <- subset(merged_data, !is.na(merged_data$PAdjust_1d_3d))
rownames(m2)=m2[,1]
m3 <- subset(merged_data, !is.na(merged_data$PAdjust_3d_5d))
rownames(m3)=m3[,1]
m4 <- subset(merged_data, !is.na(merged_data$PAdjust_5d_7d))
rownames(m4)=m4[,1]
m5 <- subset(merged_data, !is.na(merged_data$PAdjust_7d_9d))
rownames(m5)=m5[,1]

png("dge/de_genes.png", height = 17, width = 17, units = "cm", res = 300)
par(mfrow=c(2,3),mgp=c(2.2, 0.7, 0))
par(mar=c(7,7, 5.5, 3.5) - 3.0)
bl <- barplot(as.matrix(de_counts), main="A", ylab= "Number of DE transcripts", xlab= "Time intervals",ylim=c(-1200,1200), beside=TRUE, col=c("red", "blue"), names.arg=c("0d_1d","1d_3d","3d_5d","5d_7d","7d_9d"),cex.names=0.6)
text(x= bl, y= as.matrix(de_counts), labels=as.character(c(nrow(i1),nrow(d1),nrow(i2),nrow(d2),nrow(i3),nrow(d3),nrow(i4),nrow(d4),nrow(i5),nrow(d5))), xpd=TRUE, pos=c(3,1))
abline(0, 0, col = "black")

plot(m1$logCPM_0d_1d, m1$logFC_0d_1d, main="B",pch=16, ylim = c(-11,14), xlim=c(-3,18),col=rgb(0,0,0,40,maxColorValue=255), ylab = "log2FC", xlab = "log2CPM")
points(m1[m1$PAdjust_0d_1d < .05,]$logCPM_0d_1d, m1[m1$PAdjust_0d_1d <  .05,]$logFC_0d_1d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
abline(h=c(-1,1), col="grey")

plot(m2$logCPM_1d_3d, m2$logFC_1d_3d, main="C", pch=16, ylim = c(-11,14), xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), xaxt="n", ylab = "", xlab = "")
points(m2[m2$PAdjust_1d_3d < .05,]$logCPM_1d_3d, m2[m2$PAdjust_1d_3d <  .05,]$logFC_1d_3d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
abline(h=c(-1,1), col="grey")

plot(m3$logCPM_3d_5d, m3$logFC_3d_5d, main="D", pch=16, ylim = c(-11,14), xlim=c(-3,18),col=rgb(0,0,0,40,maxColorValue=255), ylab = "log2FC", xlab = "log2CPM")
points(m3[m3$PAdjust_3d_5d < .05,]$logCPM_3d_5d, m3[m3$PAdjust_3d_5d <  .05,]$logFC_3d_5d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
abline(h=c(-1,1), col="grey")

plot(m4$logCPM_5d_7d, m4$logFC_5d_7d, main="E", pch=16, ylim = c(-11,14),xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), xaxt="n", ylab = "", xlab = "")
points(m4[m4$PAdjust_5d_7d < .05,]$logCPM_5d_7d, m4[m4$PAdjust_5d_7d <  .05,]$logFC_5d_7d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
abline(h=c(-1,1), col="grey")

plot(m5$logCPM_7d_9d, m5$logFC_7d_9d, main="F", pch=16, ylim = c(-11,14),xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), xaxt="n", ylab = "", xlab = "")
points(m5[m5$PAdjust_7d_9d < .05,]$logCPM_7d_9d, m5[m5$PAdjust_7d_9d <  .05,]$logFC_7d_9d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
abline(h=c(-1,1), col="grey")
dev.off()
```

**Source:** [test_dge.r](test_dge.r)

Gene ontology
-------------

1. Parse accession numbers from annotated assembly.
2. Use [Entrez](http://www.ncbi.nlm.nih.gov/sites/gquery), [Uniprot](http://www.uniprot.org/), and [EBI
QuickGO](http://www.ebi.ac.uk/QuickGO/) to fetch list of GO ids associated with
the protein. (should I fetch InterPro?)
3. For each GO id fetch hierarchy, terms, and evidence codes.
4. Associate terms with contigs.
5. Done.
