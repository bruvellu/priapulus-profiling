#!/bin/bash

set -e

# Move to root directory.
cd ../

# Prepare data for STEM analysis.
R --save < stem_analysis/stem_prepare.r
wait
