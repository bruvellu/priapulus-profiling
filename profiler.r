profiler <- function(transcript) {
    print(transcript)

    # Create subset for transcript.
    t_subset <- subset(merged_data, reference == transcript)

    # Vector with values to be plot.
    avg_vector <- c(t_subset$avg_norm_count_0d, t_subset$avg_norm_count_1d,
                    t_subset$avg_norm_count_3d, t_subset$avg_norm_count_5d,
                    t_subset$avg_norm_count_7d, t_subset$avg_norm_count_9d)

    # Length of vector with averages.
    len_avg_vector = length(avg_vector)

    # Plot average values per time point.
    plot(avg_vector, main=transcript, ylab="Average of Normalized Counts", xlab="Time Points", pch=1, xaxt="n")
    axis(1, 1:len_avg_vector, c("0d", "1d", "3d", "5d", "7d", "9d"))

    # Segment 0d_1d.
    segments(1, avg_vector[1], 2, avg_vector[2], col=get_color(t_subset$PAdjust_0d_1d), lwd=3)
    # Segment 1d_3d.
    segments(2, avg_vector[2], 3, avg_vector[3], col=get_color(t_subset$PAdjust_1d_3d), lwd=3)
    # Segment 3d_5d.
    segments(3, avg_vector[3], 4, avg_vector[4], col=get_color(t_subset$PAdjust_3d_5d), lwd=3)
    # Segment 5d_7d.
    segments(4, avg_vector[4], 5, avg_vector[5], col=get_color(t_subset$PAdjust_5d_7d), lwd=3)
    # Segment 7d_9d.
    segments(5, avg_vector[5], 6, avg_vector[6], col=get_color(t_subset$PAdjust_7d_9d), lwd=3)
}

get_color <- function(pvalue) {
    if (pvalue < .05) {
        return("red")
    }
    else {
        return("black")
    }
}

# TODO Load profile category into merged_data.
