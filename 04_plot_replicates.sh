#!/bin/bash

set -e

# Build scatter plots for replicates.
R --save < build_scatter_plots.r
