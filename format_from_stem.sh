# Input/Output files.
INPUT=stem/genes_passing_filter
OUTPUT=stem/stem_profiles_to_r

# This script will automatically formats STEM main table of genes output. You
# can edit the file manually after running STEM using these Vim commands:
#
#   6. Remove first line: `Table of Genes Passing Filter`
#   1. Remove tab from beginning of the line with `:%s:^\t::g`
#   2. Remove initial header with `:%s:^Selected\t::g`
#   3. Remove SPOT column with `%s:\tID_\d\+\t\(\d\+\):\t\1:g`
#   4. Remove SPOT header with `:%s:SPOT\t::g`
#   5. Change count to delta in header: `:%s:avg_norm_count:avg_norm_delta:g`
#   7. Lowercase everything with: `ggVGu`
#   8. Uppercase locus and transcript on deflines with `:%s:locus_:Locus_:g` and `:%s:transcript_:Transcript_:g`

# Remove quotes.
SKIP_1ST='2,$'

# Add "transcript" to first column.
RM_TAB='s:^\t::g'

# Substitute white space for tabs.
DESELECT='s:^Selected\t::g'

# Remove SPOT column.
RM_SPOT='s:\tID_[0-9]*::g'

# Remove SPOT header.
RM_SPOTH='s:SPOT\t::g'

# Change count to delta in header.
DELTA='s:avg_norm_count:avg_norm_delta:g'

# Lowercase transcript names.
LOWER='s:LOCUS\(.*\)TRANSCRIPT:Locus\1Transcript:g'

# Execute regular expressions and write output.
sed -n "$SKIP_1ST"p $INPUT | sed "$RM_TAB" | sed "$DESELECT" | sed "$RM_SPOT" | sed "$RM_SPOTH" | sed "$DELTA" | sed "$LOWER" > $OUTPUT
