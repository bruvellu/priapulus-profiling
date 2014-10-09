multi_profiler <- function(transcripts, outdir) {
    print(transcripts)

    # Set initial maximum to 0.
    maximum <- 0
    # Count transcripts.
    n_transcripts <- length(transcripts)

    # Loop over just to get the maximum (miminum is 0).
    for (transcript in transcripts) {
        # Create subset for transcripts.
        t_subset <- subset(merged_data, reference == transcript, select=c("avg_norm_count_0d", "avg_norm_count_1d", "avg_norm_count_3d", "avg_norm_count_5d", "avg_norm_count_7d", "avg_norm_count_9d"))
        t_vector <- as.vector(t_subset, mode="numeric")
        t_max <- max(t_vector)
        if (t_max > maximum) {
            maximum <- t_max
        }
    }

    # Build file path.
    filepath <- file.path(outdir, paste0("multi_profiles", ".png"))

    # Define plot ranges.
    xrange = range(1:6)
    yrange = range(0:maximum)

    # Define plot, but only draw it with the data.
    png(filepath, width=1000, height=1000, units="px")
    plot(xrange, yrange, type="n", ylab="Average of Normalized Counts", xlab="Time Points", pch=1, xaxt="n")
    axis(1, 1:6, c("0d", "1d", "3d", "5d", "7d", "9d"))

     # Pick colors for the different lines.
    colors <- terrain.colors(n_transcripts)

    # Re-loop now to get and plot actual values.
    for (i in 1:n_transcripts) {
        t_subset <- subset(merged_data, reference == transcripts[i], select=c("avg_norm_count_0d", "avg_norm_count_1d", "avg_norm_count_3d", "avg_norm_count_5d", "avg_norm_count_7d", "avg_norm_count_9d"))
        rownames(t_subset) <- transcripts[i]
        t_vector <- as.vector(t_subset, mode="numeric")
        lines(t_vector, col=colors[i], lwd=10)
    }

    # Put a legend to identify the transcripts.
    legend("topleft", transcripts, col=colors, lwd=5)

    dev.off()

}

#plot(avg_vector, main=transcript, ylab="Average of Normalized Counts", xlab="Time Points", pch=1, xaxt="n")
#axis(1, 1:length(avg_vector), c("0d", "1d", "3d", "5d", "7d", "9d"))

## Segment 0d_1d.
#segments(1, avg_vector[1], 2, avg_vector[2], col=get_color(t_subset$PAdjust_0d_1d), lwd=3)
## Segment 1d_3d.
#segments(2, avg_vector[2], 3, avg_vector[3], col=get_color(t_subset$PAdjust_1d_3d), lwd=3)
## Segment 3d_5d.
#segments(3, avg_vector[3], 4, avg_vector[4], col=get_color(t_subset$PAdjust_3d_5d), lwd=3)
## Segment 5d_7d.
#segments(4, avg_vector[4], 5, avg_vector[5], col=get_color(t_subset$PAdjust_5d_7d), lwd=3)
## Segment 7d_9d.
#segments(5, avg_vector[5], 6, avg_vector[6], col=get_color(t_subset$PAdjust_7d_9d), lwd=3)

#text(c(1.5, 2.5, 3.5, 4.5, 5.5), c((avg_vector[1] + avg_vector[2]) / 2,
                                   #(avg_vector[2] + avg_vector[3]) / 2,
                                   #(avg_vector[3] + avg_vector[4]) / 2,
                                   #(avg_vector[4] + avg_vector[5]) / 2,
                                   #(avg_vector[5] + avg_vector[6]) / 2),
     #labels=as.character(round(log_fc, digits=1)), pos=1)

#dev.off()
