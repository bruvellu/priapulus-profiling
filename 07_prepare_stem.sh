#!/bin/bash

set -e

# Prepare data for STEM analysis.
R --save < stem_prepare.r
wait

# Format average file for STEM analysis.
sh format_to_stem.sh
wait
