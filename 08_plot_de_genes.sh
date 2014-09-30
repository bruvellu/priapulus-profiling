#!/bin/bash

set -e

# Build plots of differentially expressed genes.
R --save < build_de_plots.r
