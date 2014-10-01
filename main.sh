#!/bin/bash

# TODO Check if necessary components are installed.
# bowtie2, R, edgeR.

00_map_to_reference.sh
01_calculate_read_count.sh
02_import_to_r.sh
03_normalize.sh
04_plot_replicates.sh
05_test_expression.sh
06_plot_de_genes.sh
07_prepare_stem.sh
08_process_stem.sh
