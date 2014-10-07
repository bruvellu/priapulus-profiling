profiler <- function(transcript) {
  print(transcript)
  # Creates vector.
  avg_norm_counts <- c()
  
  # Define vector with time points.
  avg_time_points <- c("avg_norm_count_0d", "avg_norm_count_1d", "avg_norm_count_3d", "avg_norm_count_5d", "avg_norm_count_7d", "avg_norm_count_9d")
  
  for (time_point in avg_time_points) {
    print(time_point)
    value <- merged_data[merged_data$reference == transcript, time_point]
    avg_norm_counts <- append(avg_norm_counts, value)
  }
  
  return(avg_norm_counts)
}

# Build data frame with one transcript per line and column names