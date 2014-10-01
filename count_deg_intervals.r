

# Select differentially expressed genes for each interval.
count_deg <- function(df) {
    # 0d - 1d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
    i1 <- df[(df$logFC_0d_1d > 0) & (df$PAdjust_0d_1d < 0.05) & (!is.na(df$PAdjust_0d_1d)),]
    d1 <- df[(df$logFC_0d_1d < 0) & (df$PAdjust_0d_1d < 0.05) & (!is.na(df$PAdjust_0d_1d)),]

    # 1d - 3d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
    i2 <- df[(df$logFC_1d_3d > 0) & (df$PAdjust_1d_3d < 0.05) & (!is.na(df$PAdjust_1d_3d)),]
    d2 <- df[(df$logFC_1d_3d < 0) & (df$PAdjust_1d_3d < 0.05) & (!is.na(df$PAdjust_1d_3d)),]

    # 3d - 5d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
    i3 <- df[(df$logFC_3d_5d > 0) & (df$PAdjust_3d_5d < 0.05) & (!is.na(df$PAdjust_3d_5d)),]
    d3 <- df[(df$logFC_3d_5d < 0) & (df$PAdjust_3d_5d < 0.05) & (!is.na(df$PAdjust_3d_5d)),]

    # 5d - 7d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
    i4 <- df[(df$logFC_5d_7d > 0) & (df$PAdjust_5d_7d < 0.05) & (!is.na(df$PAdjust_5d_7d)),]
    d4 <- df[(df$logFC_5d_7d < 0) & (df$PAdjust_5d_7d < 0.05) & (!is.na(df$PAdjust_5d_7d)),]

    # 7d - 9d; increasing (logFC>0), significant (PAdjust<0.05), excluding NA.
    i5 <- df[(df$logFC_7d_9d > 0) & (df$PAdjust_7d_9d < 0.05) & (!is.na(df$PAdjust_7d_9d)),]
    d5 <- df[(df$logFC_7d_9d < 0) & (df$PAdjust_7d_9d < 0.05) & (!is.na(df$PAdjust_7d_9d)),]

    # Count increasing and decreasing genes per interval.
    I1 <- c(nrow(i1),-(nrow(d1)))
    I2 <- c(nrow(i2),-(nrow(d2)))
    I3 <- c(nrow(i3),-(nrow(d3)))
    I4 <- c(nrow(i4),-(nrow(d4)))
    I5 <- c(nrow(i5),-(nrow(d5)))

    # Save counts into a data frame.
    deg_counts <- data.frame(I1, I2, I3, I4, I5)

    # Return data frame.
    return(deg_counts)
}
