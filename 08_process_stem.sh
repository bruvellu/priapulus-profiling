#!/bin/bash

set -e

# Format STEM table of genes for R.
sh format_from_stem.sh
wait

# Import STEM profiles into R.
R --save < stem_import.r
wait

# Build plots of differentially expressed genes for each profile.
R --save < build_profile_plots.r
wait
