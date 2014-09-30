#!/bin/bash

set -e

# Format STEM table of genes for R.
sh format_from_stem.sh

# Import STEM profiles into R.
R --save < stem_import.r
