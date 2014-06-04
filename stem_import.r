# Put STEM data into matrix with R
stem <- read.table("stem/stem_profiles_to_r", header=TRUE)
merged_data <- merge(merged_data, stem, by.x="reference", by.y="transcript", all=TRUE)

# Save matrix in R:
write.table(merged_data, "matrix.R")
