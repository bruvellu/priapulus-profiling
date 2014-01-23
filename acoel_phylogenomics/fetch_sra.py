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

# Paired or single reads will be parsed from the field ExpXml. Examples:
# paired: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=445724
# single: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=423751

# Create list for selected IDs.
sra_selected_ids = []

# Write IDs and READ LENGTH to CSV file.
paired_sra = open('paired_sra.csv', 'w')

# Iterate through summaries and pick paired end datasets.
paired_sra.write('ID,READ_LENGTH\n')
print('\nID\tREAD_LENGTH')
for k, v in sra_summaries.iteritems():
    summary_string = v[0]['ExpXml']
    # Pattern to match: <PAIRED NOMINAL_LENGTH="200"
    re_search = re.search('<(?P<reads>PAIRED)\sNOMINAL_LENGTH="(?P<length>\d+)"', summary_string)
    try:
        reads = re_search.group('reads')
        # TODO Set minimum limit to 80 bp?
        reads_length = int(re_search.group('length'))
        sra_selected_ids.append(k)
        paired_sra.write('%s,%d\n' % (k, reads_length))
        print('%s\t%d' % (k, reads_length))
    except:
        print('%s\tdoes no satisfy parameters.' % k)
        #print('\n%s\n' % summary_string)

# Closes CSV file.
paired_sra.close()

# Just count selected IDs.
count_sra_selected_ids = len(sra_selected_ids)

print('\nTotal of %d selected IDs. %d%% of returned entries.' % (count_sra_selected_ids, count_sra_selected_ids / int(search_records['Count'])))

# TODO Use IDs to fetch the full records of these. Apparently the only return
# type is XML, so it needs to be parsed to get a tabular output, if needed.

# TODO find out how to fetch the sequences, iteratively by ID?
