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
