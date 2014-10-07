# Load plot function.
source("de_plot.r")

# Plot all data set highlighting differentially expressed genes.
#de_plots(merged_data, deg_counts_all, "plots/de_genes.png")

png("plots/de_genes.png", width=40, height=10, units="cm", res=300)
par(mfrow=c(1,6),mgp=c(3, 1, 0))
par(mar=c(3, 2, 2, 2) + 0.1)

# TODO Update de_counts for the new de_counts_all!!!
bl <- barplot(as.matrix(de_counts), main="A", ylab= "Number of DE transcripts", xlab= "Time intervals",ylim=c(-1200,1200), beside=TRUE, col=c("red", "blue"), names.arg=c("0d_1d","1d_3d","3d_5d","5d_7d","7d_9d"),cex.names=0.6)
text(x=bl, y=as.matrix(de_counts), labels=as.character(c(de_counts$I1[1],de_counts$I1[2],de_counts$I2[1],de_counts$I2[2],de_counts$I3[1],de_counts$I3[2],de_counts$I4[1],de_counts$I4[2],de_counts$I5[1],de_counts$I5[2])), xpd=TRUE, pos=c(3,1))
abline(0, 0, col = "black")

plotSmear(dge_0d1d, de.tags=merged_data[merged_data$PAdjust_0d_1d < .05,]$reference, main="B", xlim=c(-6,18), ylim=c(-11,14))
abline(h=c(-2,2), col="grey")

plotSmear(dge_1d3d, de.tags=merged_data[merged_data$PAdjust_1d_3d < .05,]$reference, main="C", xlim=c(-6,18), ylim=c(-11,14))
abline(h=c(-2,2), col="grey")

plotSmear(dge_3d5d, de.tags=merged_data[merged_data$PAdjust_3d_5d < .05,]$reference, main="D", xlim=c(-6,18), ylim=c(-11,14))
abline(h=c(-2,2), col="grey")

plotSmear(dge_5d7d, de.tags=merged_data[merged_data$PAdjust_5d_7d < .05,]$reference, main="E", xlim=c(-6,18), ylim=c(-11,14))
abline(h=c(-2,2), col="grey")

plotSmear(dge_7d9d, de.tags=merged_data[merged_data$PAdjust_7d_9d < .05,]$reference, main="F", xlim=c(-6,18), ylim=c(-11,14))
abline(h=c(-2,2), col="grey")

dev.off()