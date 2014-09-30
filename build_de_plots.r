
# Helper function to build differential expression plots.
de_plots <- function(df, filepath) {

    # Subset each interval keeping only valid p-values.
    m1 <- subset(df, !is.na(df$PAdjust_0d_1d))
    rownames(m1)=m1[,1]
    m2 <- subset(df, !is.na(df$PAdjust_1d_3d))
    rownames(m2)=m2[,1]
    m3 <- subset(df, !is.na(df$PAdjust_3d_5d))
    rownames(m3)=m3[,1]
    m4 <- subset(df, !is.na(df$PAdjust_5d_7d))
    rownames(m4)=m4[,1]
    m5 <- subset(df, !is.na(df$PAdjust_7d_9d))
    rownames(m5)=m5[,1]

    png(filepath, width=24, height=8, units="cm", res=300)
    par(mfrow=c(1,6),mgp=c(3, 1, 0))
    par(mar=c(5, 4, 4, 2) + 0.1)

    bl <- barplot(as.matrix(de_counts), main="A", ylab= "Number of DE transcripts", xlab= "Time intervals",ylim=c(-1200,1200), beside=TRUE, col=c("red", "blue"), names.arg=c("0d_1d","1d_3d","3d_5d","5d_7d","7d_9d"),cex.names=0.6)
    text(x= bl, y= as.matrix(de_counts), labels=as.character(c(nrow(i1),nrow(d1),nrow(i2),nrow(d2),nrow(i3),nrow(d3),nrow(i4),nrow(d4),nrow(i5),nrow(d5))), xpd=TRUE, pos=c(3,1))
    abline(0, 0, col = "black")

    plot(m1$logCPM_0d_1d, m1$logFC_0d_1d, main="B",pch=16, ylim=c(-11,14), xlim=c(-3,18),col=rgb(0,0,0,40,maxColorValue=255), ylab = "log2FC", xlab = "log2CPM")
    points(m1[m1$PAdjust_0d_1d < .05,]$logCPM_0d_1d, m1[m1$PAdjust_0d_1d <  .05,]$logFC_0d_1d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    plot(m2$logCPM_1d_3d, m2$logFC_1d_3d, main="C", pch=16, ylim=c(-11,14), xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), xaxt="n", ylab = "", xlab = "")
    points(m2[m2$PAdjust_1d_3d < .05,]$logCPM_1d_3d, m2[m2$PAdjust_1d_3d <  .05,]$logFC_1d_3d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    plot(m3$logCPM_3d_5d, m3$logFC_3d_5d, main="D", pch=16, ylim=c(-11,14), xlim=c(-3,18),col=rgb(0,0,0,40,maxColorValue=255), ylab = "log2FC", xlab = "log2CPM")
    points(m3[m3$PAdjust_3d_5d < .05,]$logCPM_3d_5d, m3[m3$PAdjust_3d_5d <  .05,]$logFC_3d_5d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    plot(m4$logCPM_5d_7d, m4$logFC_5d_7d, main="E", pch=16, ylim=c(-11,14),xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), xaxt="n", ylab = "", xlab = "")
    points(m4[m4$PAdjust_5d_7d < .05,]$logCPM_5d_7d, m4[m4$PAdjust_5d_7d <  .05,]$logFC_5d_7d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    plot(m5$logCPM_7d_9d, m5$logFC_7d_9d, main="F", pch=16, ylim=c(-11,14),xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), xaxt="n", ylab = "", xlab = "")
    points(m5[m5$PAdjust_7d_9d < .05,]$logCPM_7d_9d, m5[m5$PAdjust_7d_9d <  .05,]$logFC_7d_9d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    dev.off()
}

# Plots with all data.
de_plots(merged_data, "plots/de_genes.png")

# Plot Profile 8 genes.
profile8 <- read.table("stem/profile_8_to_r", header=TRUE)
stem8 <- merged_data[merged_data$reference %in% profile8$transcript,]
de_plots(stem8, "plots/de_genes_stem_8.png")

# Plot Profile 39
profile39 <- read.table("stem/profile_39_to_r", header=TRUE)
stem39 <- merged_data[merged_data$reference %in% profile39$transcript,]
de_plots(stem39, "plots/de_genes_stem_39.png")

# Plot Profile 22
profile22 <- read.table("stem/profile_22_to_r", header=TRUE)
stem22 <- merged_data[merged_data$reference %in% profile22$transcript,]
de_plots(stem22, "plots/de_genes_stem_22.png")

# Plot Profile 31 genes.
profile31 <- read.table("stem/profile_31_to_r", header=TRUE)
stem31 <- merged_data[merged_data$reference %in% profile31$transcript,]
de_plots(stem31, "plots/de_genes_stem_31.png")

# Plot Profile 17 genes.
profile17 <- read.table("stem/profile_17_to_r", header=TRUE)
stem17 <- merged_data[merged_data$reference %in% profile17$transcript,]
de_plots(stem17, "plots/de_genes_stem_17.png")

# Plot Profile 18 genes.
profile18 <- read.table("stem/profile_18_to_r", header=TRUE)
stem18 <- merged_data[merged_data$reference %in% profile18$transcript,]
de_plots(stem18, "plots/de_genes_stem_18.png")

# Plot Profile 1 genes.
profile1 <- read.table("stem/profile_1_to_r", header=TRUE)
stem1 <- merged_data[merged_data$reference %in% profile1$transcript,]
de_plots(stem1, "plots/de_genes_stem_1.png")

# Plot Profile 25 genes.
profile25 <- read.table("stem/profile_25_to_r", header=TRUE)
stem25 <- merged_data[merged_data$reference %in% profile25$transcript,]
de_plots(stem25, "plots/de_genes_stem_25.png")

#profile_6_to_r
#profile_7_to_r
#profile_10_to_r
#profile_15_to_r
