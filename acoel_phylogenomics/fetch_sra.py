#!/usr/bin/env python
'''Search & Fetch records from NCBI\'s Sequence Read Archive.'''

import argparse
import re

from Bio import Entrez


class SRADatabase:
    '''General information about SRA Database'''

    def __init__(self):
        einfo_handle = Entrez.einfo(db='sra')
        einfo = Entrez.read(einfo_handle, validate=False)

        # Define attributes.
        self.count = einfo['DbInfo']['Count']
        self.last_update = einfo['DbInfo']['LastUpdate']
        self.menu_name = einfo['DbInfo']['MenuName']
        self.description = einfo['DbInfo']['Description']
        self.link_list = einfo['DbInfo']['LinkList']
        self.field_list = einfo['DbInfo']['FieldList']
        self.einfo = einfo

        #print(self.count, self.last_update, self.menu_name, self.description,
              #self.link_list, self.field_list, self.einfo)

    # TODO Use def __unicode__ or __str__ to identify class objects.

class SRASearch:
    '''Perform search and keep IDs of SRA packages.

    Example of query:

        ((((("strategy rna seq"[Properties]) AND "platform illumina"[Properties])
        AND metazoa[Organism]) NOT vertebrata[Organism]) NOT insects[Organism]) AND
        ("2000/01/01"[Modification Date] : "3000"[Modification Date])
    '''

    def __init__(self, query, retmax, email):
        # Required arguments.
        self.query = query
        if int(retmax) > 100000:
            # Limit defined by Entrez.
            self.retmax = 100000
        else:
            self.retmax = retmax
        Entrez.email = email

        # Search metadata.
        self.count = None
        self.retstart = None
        self.query_translation = None
        self.idlist = None

        # Additional attributes.
        self.results = None
        self.database = SRADatabase()
        # TODO Add timestamp.

    def esearch(self):
        '''Search SRA packages with Entrez using query.'''
        handle = Entrez.esearch(db='sra', term=self.query, retmax=self.retmax)
        self.results = Entrez.read(handle)
        self.parse_results()
        return self.results

    def parse_results(self):
        '''Populate class attributes by parsing results.'''
        self.count = self.results['Count']
        self.retstart = self.results['RetStart']
        self.query_translation = self.results['QueryTranslation']
        self.idlist = self.results['IdList']

        #print(self.count, self.retstart, self.query_translation, self.idlist)

class SRAPackage:
    '''Fetch and store metadata from a SRA package.'''

    def __init__(self, sra_id):
        self.id = sra_id
        self.record = None

        self.accession = None
        self.library_strategy = None
        self.library_layout = None
        self.instrument_model = None
        self.taxon_id = None
        self.scientific_name = None
        self.runs = None
        self.nreads = None
        self.average = None

        self.ncbi_taxid = None

        self.efetch()

    def efetch(self):
        '''Fetch package metadata from Entrez'''
        handle = Entrez.efetch(db='sra', id=self.id)
        self.record = handle.read()
        self.extract()

    def extract(self):
        '''Extract relevant fields from summary.'''

        # Fields with attributes.
        fields = {}

        # Fields to be parsed.
        regexes = {
            'accession': '<EXPERIMENT\s+.*?accession="(?P<accession>.*?)".*?>',
            'library_strategy': '<LIBRARY_STRATEGY>(?P<library_strategy>.*?)<\/LIBRARY_STRATEGY>',
            'library_layout': '<LIBRARY_LAYOUT>\s*<(?P<type>SINGLE|PAIRED)',
            'instrument_model': '<INSTRUMENT_MODEL>(?P<instrument_model>.*?)<\/INSTRUMENT_MODEL>',
            'taxon_id': '<TAXON_ID>(?P<taxon_id>.*?)<\/TAXON_ID>',
            'scientific_name': '<SCIENTIFIC_NAME>(?P<scientific_name>.*?)<\/SCIENTIFIC_NAME>',
            'runs': '<RUN\s+.*?accession="(?P<accession>.*?)"\s+.*?total_spots="(?P<total_spots>.*?)"\s+.*?total_bases="(?P<total_bases>.*?)"\s+.*?size="(?P<size>.*?)"\s+.*?published="(?P<published>.*?)"\s+.*?>',
            'nreads': '<Statistics\s+.*?nreads="(?P<nreads>.*?)"\s+.*?>',
            'average': '<Read\s+.*?average="(?P<average>.*?)"\s+.*?\/>',
        }

        # Iterate over regexes to parse attributes.
        # TODO handle multiple matches like "runs", "nreads", and "average".
        for field, regex in regexes.iteritems():
            re_search = re.search(regex, self.record)
            if re_search and re_search.groups(0):
                fields[field] = re_search.groups(0)[0]
            else:
                fields[field] = ''

        print('')
        print(fields)


class FilterPackages:
    '''Filter results based on package metadata.'''
    # TODO Plan a way to effectively filter results. Maybe get scipy or pandas
    # help?

    def __init__(self, packages):
        self.packages = packages


def main():
    '''Parse arguments and call SRA search.

    Main function simply parses arguments from command line input and assures
    everything is ok to instantiate the SRA search class.
    '''

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Search & Fetch records from NCBI\'s Sequence Read Archive.',
                                     epilog='Work out those reads, dude.')
    parser.add_argument('-s', '--search',
                        help='put search terms between "quotes"',
                        type=str, required=True)
    parser.add_argument('-m', '--maximum',
                        help='maximum number of records to be retrieved',
                        default='20')
    parser.add_argument('-o', '--output',
                        help='indicate output CSV file',
                        required=True)
    parser.add_argument('-e', '--email',
                        help='an email address is required for Entrez',
                        required=True)
    args = parser.parse_args()

    # Instantiate search object.
    sra_search = SRASearch(query=args.search, retmax=args.maximum,
                           email=args.email)

    # Execute search itself.
    sra_search.esearch()

    # Fetch metadata from packages.
    packages = [SRAPackage(sra_id) for sra_id in sra_search.idlist]

if __name__ == '__main__':
    main()


## Paired or single reads will be parsed from the field ExpXml. Examples:
## paired: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=445724
## single: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=423751

## Create list for selected IDs.
#sra_selected_ids = []

## Write IDs and READ LENGTH to CSV file.
#paired_sra = open('paired_sra.csv', 'w')

## Iterate through summaries and pick paired end datasets.
#header = ['SRA_ID', 'ORGANISM', 'NCBI_TAXID', 'STRATEGY', 'PLATFORM', 'TYPE',
          #'READ_LENGTH', 'READS', 'SPOTS', 'BASES']

#paired_sra.write(', '.join(header) + '\n')
#print('\n' + '\t '.join(header))
#for k, v in sra_summaries.iteritems():

    ## Pattern to match: <PAIRED NOMINAL_LENGTH="200"
    #summary_string = v[0]['ExpXml']

    ## Extract fields from record summary.
    #fields = extract_summary(summary_string)

    #sra_selected_ids.append(k)
    #row = [k, fields['organism'], fields['ncbi_taxid'], fields['strategy'],
            #fields['platform'], fields['type'] ,fields['length'],
            #fields['reads'], fields['spots'], fields['spots']]
    #if fields['type'] == 'PAIRED':
        #paired_sra.write(', '.join(row) + '\n')
    #print('\t '.join(row))

## Closes CSV file.
#paired_sra.close()

## Just count selected IDs.
#count_sra_selected_ids = len(sra_selected_ids)

#print('\nTotal of %d selected IDs. %d%% of returned entries.' % (count_sra_selected_ids, count_sra_selected_ids / int(search_records['Count'])))
