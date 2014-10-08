# Load edgeR
library(edgeR)

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

# Calculate the average between replicates.
avg_norm_count_0d <- (merged_data$norm_count_Pc1_0d + merged_data$norm_count_Pc2_0d) / 2
avg_norm_count_1d <- (merged_data$norm_count_Pc1_1d + merged_data$norm_count_Pc2_1d) / 2
avg_norm_count_3d <- (merged_data$norm_count_Pc1_3d + merged_data$norm_count_Pc2_3d) / 2
avg_norm_count_5d <- (merged_data$norm_count_Pc1_5d + merged_data$norm_count_Pc2_5d) / 2
avg_norm_count_7d <- merged_data$norm_count_Pc2_7d
avg_norm_count_9d <- (merged_data$norm_count_Pc1_9d + merged_data$norm_count_Pc2_9d) / 2

# Build data frame for average values.
avg_data_norm <- cbind(avg_norm_count_0d, avg_norm_count_1d, avg_norm_count_3d, avg_norm_count_5d, avg_norm_count_7d, avg_norm_count_9d)

# Rename columns.
colnames(avg_data_norm) <- c("avg_norm_count_0d", "avg_norm_count_1d", "avg_norm_count_3d", "avg_norm_count_5d", "avg_norm_count_7d", "avg_norm_count_9d")

# Bind average counts with matrix.
merged_data <- cbind(merged_data, avg_data_norm)
