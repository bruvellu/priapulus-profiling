# Build and save scatter plots for replicates.

pdf("plots/Pc1_0d__Pc2_0d.pdf")
plot(merged_data$norm_count_Pc1_0d, merged_data$norm_count_Pc2_0d, log="xy", xlab="Pc1_0d normalized counts", ylab="Pc2_0d normalized counts", main="0d log distribution of normalized counts")
dev.off()

pdf("plots/Pc1_1d__Pc2_1d.pdf")
plot(merged_data$norm_count_Pc1_1d, merged_data$norm_count_Pc2_1d, log="xy", xlab="Pc1_1d normalized counts", ylab="Pc2_1d normalized counts", main="1d log distribution of normalized counts")
dev.off()

pdf("plots/Pc1_3d__Pc2_3d.pdf")
plot(merged_data$norm_count_Pc1_3d, merged_data$norm_count_Pc2_3d, log="xy", xlab="Pc1_3d normalized counts", ylab="Pc2_3d normalized counts", main="3d log distribution of normalized counts")
dev.off()

pdf("plots/Pc1_5d__Pc2_5d.pdf")
plot(merged_data$norm_count_Pc1_5d, merged_data$norm_count_Pc2_5d, log="xy", xlab="Pc1_5d normalized counts", ylab="Pc2_5d normalized counts", main="5d log distribution of normalized counts")
dev.off()

#pdf("plots/Pc1_7d__Pc2_7d.pdf")
#plot(merged_data$norm_count_Pc1_7d, merged_data$norm_count_Pc2_7d, log="xy", xlab="Pc1_7d normalized counts", ylab="Pc2_7d normalized counts", main="7d log distribution of normalized counts")
#dev.off()

pdf("plots/Pc1_9d__Pc2_9d.pdf")
plot(merged_data$norm_count_Pc1_9d, merged_data$norm_count_Pc2_9d, log="xy", xlab="Pc1_9d normalized counts", ylab="Pc2_9d normalized counts", main="9d log distribution of normalized counts")
dev.off()

png("plots/rep_plots.png", width=1250, height=250)
par(mfrow=c(1,5))
plot(merged_data$norm_count_Pc1_0d, merged_data$norm_count_Pc2_0d, log="xy", asp=1, xlab="Pc1_0d normalized counts", ylab="Pc2_0d normalized counts", main="0d log distribution of normalized counts")
plot(merged_data$norm_count_Pc1_1d, merged_data$norm_count_Pc2_1d, log="xy", asp=1, xlab="Pc1_1d normalized counts", ylab="Pc2_1d normalized counts", main="1d log distribution of normalized counts")
plot(merged_data$norm_count_Pc1_3d, merged_data$norm_count_Pc2_3d, log="xy", asp=1, xlab="Pc1_3d normalized counts", ylab="Pc2_3d normalized counts", main="3d log distribution of normalized counts")
plot(merged_data$norm_count_Pc1_5d, merged_data$norm_count_Pc2_5d, log="xy", asp=1, xlab="Pc1_5d normalized counts", ylab="Pc2_5d normalized counts", main="5d log distribution of normalized counts")
plot(merged_data$norm_count_Pc1_9d, merged_data$norm_count_Pc2_9d, log="xy", asp=1, xlab="Pc1_9d normalized counts", ylab="Pc2_9d normalized counts", main="9d log distribution of normalized counts")
dev.off()
