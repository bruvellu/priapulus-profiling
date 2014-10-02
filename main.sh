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

time sh 00_map_to_reference.sh
time sh 01_calculate_read_count.sh
time sh 02_import_to_r.sh
time sh 03_normalize.sh
time sh 04_plot_replicates.sh
time sh 05_test_expression.sh
time sh 06_plot_de_genes.sh
