# Put STEM data into matrix with R
stem <- read.table("stem/stem_profiles_to_r", header=TRUE)
merged_data <- merge(merged_data, stem, by.x="reference", by.y="transcript", all=TRUE)

# Save matrix in R:
write.table(merged_data, "matrix.R")

# Load function to count.
source("count_deg_intervals.r")

# Counts for stem8.
profile8 <- read.table("stem/profile_8_to_r", header=TRUE)
stem8 <- merged_data[merged_data$reference %in% profile8$transcript,]
deg_counts_stem8 <- count_deg(stem8)

# Counts for stem39.
profile39 <- read.table("stem/profile_39_to_r", header=TRUE)
stem39 <- merged_data[merged_data$reference %in% profile39$transcript,]
deg_counts_stem39 <- count_deg(stem39)

# Counts for stem22.
profile22 <- read.table("stem/profile_22_to_r", header=TRUE)
stem22 <- merged_data[merged_data$reference %in% profile22$transcript,]
deg_counts_stem22 <- count_deg(stem22)

# Counts for stem31.
profile31 <- read.table("stem/profile_31_to_r", header=TRUE)
stem31 <- merged_data[merged_data$reference %in% profile31$transcript,]
deg_counts_stem31 <- count_deg(stem31)

# Counts for stem17.
profile17 <- read.table("stem/profile_17_to_r", header=TRUE)
stem17 <- merged_data[merged_data$reference %in% profile17$transcript,]
deg_counts_stem17 <- count_deg(stem17)

# Counts for stem18.
profile18 <- read.table("stem/profile_18_to_r", header=TRUE)
stem18 <- merged_data[merged_data$reference %in% profile18$transcript,]
deg_counts_stem18 <- count_deg(stem18)

# Counts for stem1.
profile1 <- read.table("stem/profile_1_to_r", header=TRUE)
stem1 <- merged_data[merged_data$reference %in% profile1$transcript,]
deg_counts_stem1 <- count_deg(stem1)

# Counts for stem25.
profile25 <- read.table("stem/profile_25_to_r", header=TRUE)
stem25 <- merged_data[merged_data$reference %in% profile25$transcript,]
deg_counts_stem25 <- count_deg(stem25)

