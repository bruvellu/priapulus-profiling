# Prepare data for STEM analysis.
R --save < stem_prepare.r

# Format average file for STEM analysis.
sh format_to_stem.sh

# Wait?

# Format STEM table of genes for R.
sh format_from_stem.sh

# Import STEM profiles into R.
R --save < stem_import.r
