# Input/Output files.
INPUT=$1
OUTPUT=$1_to_r

# Remove quotes.
SKIP='5,$'

# Substitute white space for tabs.
DESELECT='s:^Selected.*$:transcript:g'

# Lowercase transcript names.
LOWER='s:^.*LOCUS\(.*\)TRANSCRIPT_\([0-9]*\).*$:Locus\1Transcript_\2:g'

# Execute regular expressions and write output.
sed -n "$SKIP"p $INPUT | sed "$DESELECT" | sed "$LOWER" > $OUTPUT
