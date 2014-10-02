#!/bin/bash

set -e

# Move to root directory.
cd ../

# Import STEM profiles into R.
R --save < stem_analysis/stem_import.r
wait

# Build plots of differentially expressed genes for each profile.
R --save < stem_analysis/build_profile_plots.r
wait
