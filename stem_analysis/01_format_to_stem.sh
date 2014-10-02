#!/bin/bash

set -e

# Format average file for STEM analysis.
# Input/Output files.
INPUT=avg_data_norm.tsv
OUTPUT=avg_stem_input

# This script will automatically formats average counts for STEM analysis. Note
# that you can also edit the file manually before running STEM using these Vim
# commands:
#
#	1. Manually add `'"transcript" '` as the first column name (without single quotes, note the white space).
#	2. Remove quotes with `:%s:"::g`.
#	3. Substitute white space for tabs `:%s:\s\+:\t:g`.

# Remove quotes.
RM_QUOTES='s:"::g'

# Add "transcript" to first column.
TRANSCRIPT='s:avg_norm_count_0d:transcript avg_norm_count_0d:'

# Substitute white space for tabs.
TABS='s:\s+:\t:g'

# Execute regular expressions and write output.
sed "$RM_QUOTES" $INPUT | sed "$TRANSCRIPT" | sed -r "$TABS" > $OUTPUT
