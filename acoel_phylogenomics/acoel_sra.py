#!/usr/bin/env python
'''Search, fetch, and filter SRA packages for acoels phylogenomics.'''

from fetch_sra import SRASearch, SRAPackage, FilterPackages
from email import email_bruno

# Define basename for output files.
basename = 'sra_paired_gte70bp'
#basename = 'testing'

# Search for RNAseq by Illumina or 454 only metazoans, but not vertebrates nor insects
# later than year 2000.
query = '''((((((strategy rna seq[Properties]) AND platform
illumina[Properties]) OR platform LS454[Properties]) AND metazoa[Organism]) NOT
vertebrata[Organism]) NOT insects[Organism]) AND ("200"[Modification Date] : "3000"[Modification Date])'''

# Maximum number of returned results.
retmax = 5000

# Instantiate search object.
sra_search = SRASearch(query=query, retmax=retmax, email=email_bruno)

# Execute search itself.
sra_search.esearch()

# Fetch metadata from packages.
packages = [SRAPackage(sra_id) for sra_id in sra_search.idlist]

# Store packages in data frame for filtering.
package_filter = FilterPackages(packages)

# copy working data frame.
df = package_filter.data_frame

# Filter booleans.
filtered_df = df[df.library_layout == 'PAIRED'][df.nreads > 1][df.read_average >= 70]

# Sort data buy lineage.
sorted_df = filtered_df.sort('lineage')

# Write CSV out.
package_filter.filtered_data_frame = sorted_df
package_filter.write_csv(basename + '.csv')

# Write unique list of taxa.
unique = package_filter.filtered_data_frame.lineage.unique()
unique.tofile('unique_' + basename + '.txt', sep='\n')
