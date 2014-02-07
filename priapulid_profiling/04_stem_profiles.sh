# Prepare data for STEM analysis.
R --save < stem_prepare.r

# Note that you have to manually edit the file `average` before running STEM. Vim commands:
#
#	1. Manually add `'"transcript" '` as the first column name (without single quotes, note the white space).
#	2. Remove quotes with `:%s:"::g`.
#	3. Substitute white space for tabs `:%s:\s\+:\t:g`.

# Import STEM profiles into R.
R --save < stem_import.r
