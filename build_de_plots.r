# Load plot function.
source("de_plot.r")

# Plot all data set highlighting differentially expressed genes.
de_plots(merged_data, deg_counts_all, "plots/de_genes.png")
