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
