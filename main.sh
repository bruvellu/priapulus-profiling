#!/bin/bash

# TODO Check if necessary components are installed.
# bowtie2, R, edgeR.

# edgeR
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#Using Bioconductor version 2.14 (BiocInstaller 1.14.2), R version
#  3.1.1.
#trying URL 'http://bioconductor.org/packages/2.14/bioc/src/contrib/limma_3.20.9.tar.gz'
#trying URL 'http://bioconductor.org/packages/2.14/bioc/src/contrib/edgeR_3.6.8.tar.gz'

sh 00_map_to_reference.sh
sh 01_calculate_read_count.sh
sh 02_import_to_r.sh
sh 03_normalize.sh
sh 04_plot_replicates.sh
sh 05_test_expression.sh
sh 06_plot_de_genes.sh
sh 07_prepare_stem.sh
sh 08_process_stem.sh
