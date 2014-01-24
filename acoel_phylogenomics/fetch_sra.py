import re
from Bio import Entrez

'''Utility to filter and download SRA records.'''

def extract_summary(summary):
    '''Use regex to parse summary fields.
    <Summary>
        <Title>Transcriptome Analysis of house mouse</Title>
        <Platform instrument_model="Illumina HiSeq 2000">ILLUMINA</Platform>
        <Statistics total_runs="1" total_spots="12891427" total_bases="1933714050" total_size="1062140634" load_done="true" cluster_name="public"/>
    </Summary>

    <Submitter acc="ERA237533" center_name="SC" contact_name="" lab_name=""/>
    <Experiment acc="ERX284729" ver="1" status="public" name=""/>
    <Study acc="ERP002100" name="Global_gene_expression_profiles_during_differentiation_of_mouse_embryonic_stem_cells_to_macrophages_"/>
    <Organism taxid="10090" CommonName="mus musculus"/>
    <Sample acc="ERS222189" name="mus musculus"/>
    <Instrument ILLUMINA="Illumina HiSeq 2000"/>

    <Library_descriptor>
        <LIBRARY_NAME>6884103</LIBRARY_NAME>
        <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
            <PAIRED NOMINAL_LENGTH="300"/>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Illumina cDNA protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
    </Library_descriptor>

    <Biosample id="2229746" acc="SAMEA2058346" sample_id="556692" sample_acc="ERS222189"/>
    <Bioproject>ERP002100</Bioproject>
    '''

    # Store fields in dictionary.
    summary_fields = {}

    # Fields to be parsed.
    fields_regexes = {
        'type'       : '<LIBRARY_LAYOUT>\s*<(?P<type>SINGLE|PAIRED)',
        'length'     : 'NOMINAL_LENGTH="(?P<length>\d+)"',
        'strategy'   : '<LIBRARY_STRATEGY>(?P<strategy>.*)</LIBRARY_STRATEGY>',
        'platform'   : 'instrument_model="(?P<platform>[\w\s]*)"',
        'spots'      : 'total_spots="(?P<spots>\d+)"',
        'bases'      : 'total_bases="(?P<bases>\d+)"',
        'reads'      : 'total_size="(?P<reads>\d+)"',
        'ncbi_taxid' : 'taxid="(?P<ncbi_taxid>\d+)"',
        'organism'   : 'ScientificName="(?P<organism>[\w\s]*)"',
    }

    for field, regex in fields_regexes.iteritems():
        re_search = re.search(regex, summary_string)
        if re_search and re_search.groups(0):
            #print(re_search.groups(0))
            summary_fields[field] = re_search.groups(0)[0]
        else:
            summary_fields[field] = ''

    #print(summary_fields)

    return summary_fields


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
header = ['SRA_ID', 'ORGANISM', 'NCBI_TAXID', 'STRATEGY', 'PLATFORM', 'TYPE',
          'READ_LENGTH', 'READS', 'SPOTS', 'BASES']

paired_sra.write(', '.join(header) + '\n')
print('\n' + '\t '.join(header))
for k, v in sra_summaries.iteritems():

    # Pattern to match: <PAIRED NOMINAL_LENGTH="200"
    summary_string = v[0]['ExpXml']

    # Extract fields from record summary.
    fields = extract_summary(summary_string)

    sra_selected_ids.append(k)
    row = [k, fields['organism'], fields['ncbi_taxid'], fields['strategy'],
            fields['platform'], fields['type'] ,fields['length'],
            fields['reads'], fields['spots'], fields['spots']]
    if fields['type'] == 'PAIRED':
        paired_sra.write(', '.join(row) + '\n')
    print('\t '.join(row))

# Closes CSV file.
paired_sra.close()

# Just count selected IDs.
count_sra_selected_ids = len(sra_selected_ids)

print('\nTotal of %d selected IDs. %d%% of returned entries.' % (count_sra_selected_ids, count_sra_selected_ids / int(search_records['Count'])))

# TODO Use IDs to fetch the full records of these. Apparently the only return
# type is XML, so it needs to be parsed to get a tabular output, if needed.

# TODO find out how to fetch the sequences, iteratively by ID?
