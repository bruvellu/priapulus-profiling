#!/usr/bin/env python
'''Utility to filter and download SRA records.'''

import argparse
import re
import sys

from Bio import Entrez


# TODO Adapt this function to full record inside SRAPackage class.
def extract_summary(summary_string):
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
        self.summary = None
        self.fields = None

        self.efetch()

    def efetch(self):
        '''Fetch package metadata from Entrez'''
        # TODO Use efetch instead of esummary to get read_length correctly.
        handle = Entrez.esummary(db='sra', id=self.id)
        self.summary = Entrez.read(handle)
        self.extract()

    def extract(self):
        '''Extract relevant fields from summary.'''
        # TODO Make a real function and define attribute fields individually
        # from full record metadata.
        string = self.summary[0]['ExpXml']
        self.fields = extract_summary(string)

        #print(self.fields)


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
