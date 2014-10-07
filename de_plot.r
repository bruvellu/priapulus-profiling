# Helper function to build differential expression plots.
de_plots <- function(df, de_counts, filepath) {

    png(filepath, width=40, height=10, units="cm", res=300)
    par(mfrow=c(1,6),mgp=c(3, 1, 0))
    par(mar=c(3, 2, 2, 2) + 0.1)

    bl <- barplot(as.matrix(de_counts), main="A", ylab= "Number of DE transcripts", xlab= "Time intervals",ylim=c(-1200,1200), beside=TRUE, col=c("red", "blue"), names.arg=c("0d_1d","1d_3d","3d_5d","5d_7d","7d_9d"),cex.names=0.6)
    text(x=bl, y=as.matrix(de_counts), labels=as.character(c(de_counts$I1[1],de_counts$I1[2],de_counts$I2[1],de_counts$I2[2],de_counts$I3[1],de_counts$I3[2],de_counts$I4[1],de_counts$I4[2],de_counts$I5[1],de_counts$I5[2])), xpd=TRUE, pos=c(3,1))
    abline(0, 0, col = "black")

    plot(df$logCPM_0d_1d, df$logFC_0d_1d, main="B",pch=16, ylim=c(-11,14), xlim=c(-3,18),col=rgb(0,0,0,40,maxColorValue=255), ylab = "log2FC", xlab = "log2CPM")
    points(df[df$PAdjust_0d_1d < .05,]$logCPM_0d_1d, df[df$PAdjust_0d_1d <  .05,]$logFC_0d_1d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    plot(df$logCPM_1d_3d, df$logFC_1d_3d, main="C", pch=16, ylim=c(-11,14), xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), ylab = "", xlab = "")
    points(df[df$PAdjust_1d_3d < .05,]$logCPM_1d_3d, df[df$PAdjust_1d_3d <  .05,]$logFC_1d_3d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    plot(df$logCPM_3d_5d, df$logFC_3d_5d, main="D", pch=16, ylim=c(-11,14), xlim=c(-3,18),col=rgb(0,0,0,40,maxColorValue=255), ylab = "log2FC", xlab = "log2CPM")
    points(df[df$PAdjust_3d_5d < .05,]$logCPM_3d_5d, df[df$PAdjust_3d_5d <  .05,]$logFC_3d_5d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    plot(df$logCPM_5d_7d, df$logFC_5d_7d, main="E", pch=16, ylim=c(-11,14),xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), ylab = "", xlab = "")
    points(df[df$PAdjust_5d_7d < .05,]$logCPM_5d_7d, df[df$PAdjust_5d_7d <  .05,]$logFC_5d_7d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    plot(df$logCPM_7d_9d, df$logFC_7d_9d, main="F", pch=16, ylim=c(-11,14),xlim=c(-3,18), col=rgb(0,0,0,40,maxColorValue=255), ylab = "", xlab = "")
    points(df[df$PAdjust_7d_9d < .05,]$logCPM_7d_9d, df[df$PAdjust_7d_9d <  .05,]$logFC_7d_9d, pch=16, col=rgb(255,0,0,70,maxColorValue=255))
    abline(h=c(-1,1), col="grey")

    dev.off()
}
