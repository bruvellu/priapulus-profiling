#!/bin/bash

# TODO Check if necessary components are installed.
# bowtie2, R, edgeR.

00_map_to_reference.sh
01_calculate_read_count.sh
02_import_to_r.sh
03_normalize.sh
04_plot_replicates.sh
05_stem_profiles.sh
06_import_stem.sh
07_test_expression.sh
08_plot_de_genes.sh
