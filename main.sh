#!/bin/bash

# TODO Check if necessary components are installed.
# bowtie2, R, edgeR.

sh 00_map_to_reference.sh
sh 01_calculate_read_count.sh
sh 02_import_to_r.sh
sh 03_normalize.sh
sh 04_plot_replicates.sh
sh 05_test_expression.sh
sh 06_plot_de_genes.sh
sh 07_prepare_stem.sh
sh 08_process_stem.sh
