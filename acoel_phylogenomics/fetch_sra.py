import re
from Bio import Entrez

'''Utility to filter and download SRA records.'''

# Needed email.
Entrez.email = 'organelas@gmail.com'

# Terms for query.
search_terms = '((((("strategy rna seq"[Properties]) AND "platform illumina"[Properties]) AND metazoa[Organism]) NOT vertebrata[Organism]) NOT insects[Organism]) AND ("2000/01/01"[Modification Date] : "3000"[Modification Date])'

# First search for records matching the terms.
search_handle = Entrez.esearch(db='sra', term=search_terms, retmax=20)
search_records = Entrez.read(search_handle)

# Print some information for the records.
print('Found %s entries' % search_records['Count'])
print('Search terms: %s' % search_records['QueryTranslation'])

# Save list with IDs.
sra_ids = search_records['IdList']

# Dictionary for summaries.
sra_summaries = {}

# Fetch basic infos to filter entries.
for sra_id in sra_ids:
    summary_handle = Entrez.esummary(db='sra', id=sra_id)
    summary_record = Entrez.read(summary_handle)
    sra_summaries[sra_id] = summary_record

# Create list for selected IDs.
sra_selected_ids = []

# Iterate through summaries and pick paired end datasets.
print('\nID\tREAD LENGTH')
for k, v in sra_summaries.iteritems():
    summary_string = v[0]['ExpXml']
    # Pattern to match: <PAIRED NOMINAL_LENGTH="200"
    re_search = re.search('<(?P<reads>PAIRED)\sNOMINAL_LENGTH="(?P<length>\d+)"', summary_string)
    try:
        reads = re_search.group('reads')
        # TODO Set minimum limit to 80 bp?
        reads_length = int(re_search.group('length'))
        print('%s\t%d' % (k, reads_length))
        sra_selected_ids.append(k)
    except:
        print('%s\tdoes no satisfy parameters.' % k)
        print('\n%s\n' % summary_string)

print('\nTotal of %d selected IDs.' % len(sra_selected_ids))

# TODO Parse read information, single or paired end and read length.
# Probably just regex the ExpXml string, Biopython does not parse this.
# Examples for:
# paired: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=445724
# single: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=423751

# TODO Select only archives with paired reads.

# TODO Use IDs to fetch the full records of these. Apparently the only return
# type is XML, so it needs to be parsed to get a tabular output, if needed.

# TODO find out how to fetch the sequences, iteratively by ID?

