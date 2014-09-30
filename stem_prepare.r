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
write.table(avg_data_norm, "stem/avg_data_norm.tsv")
