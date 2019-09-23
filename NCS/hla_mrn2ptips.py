#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 08:50:58 2019
Updated on

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

hla_mrn2ptips.py
    Reads the Lab*.csv file of HLA information and 
    replaces the MRNs with Phenotips IDs

v 1.0 - initial code version.
    
"""
__author__ = 'Susan Huse'
__date__ = 'June 27, 2019'
__version__ = '1.0'
__copyright__ = 'No copyright protection, can be used freely'

import sys
import os
import re
import datetime
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from ncbr_huse import send_update, err_out, test_file
#import json
from ncbr_bsi import read_conf, send_curl, get_bsi_session, get_bsi_name, bsi_query

######################################
#   
# Output the dbGaP Dataset (DS) files
#
######################################

def filter_results(mrnDF, batch):
    # only samples up through the specified batch number
    mrnDF = mrnDF[mrnDF['Batch_Received'] <= batch]

    # only data that are
    mrnDF = mrnDF[mrnDF['Batch_Received'].str.contains("BATCH")]

    # get rid of canceled orders
    mrnDF = mrnDF.loc[~mrnDF['CRIS_Order_Status'].str.contains('Cancel')]
    mrnDF = mrnDF[mrnDF['CRIS_Order_Status'] != 'Auto Complete']

    return(mrnDF)
    
######################################
#   
# Output the dbGaP Dataset (DS) files
#
######################################

def write_files(mrnDF, filenames):
    #filenames 0: file_consent, 1:file_mapping, 2:file_ped, 3:file_pheno, 4:file_samples):

    # Consent - ID and Consent = 1
    ##SUE!! Need to get the actual new dbGaP IDs
    mrnDF['SUBJECT_ID'] = 'New_dbGap_ID'
    mrnDF['CONSENT'] = 1
    mrnDF['SUBJECT_SOURCE'] = 'NIAID CSI'
    mrnDF['SOURCE_SUBJECT_ID'] = mrnDF['Phenotips_ID']
    #consentFields = ['SUBJECT_ID', 'CONSENT', 'SUBJECT_SOURCE', 'SOURCE_SUBJECT_ID']
    consentFields = ['SUBJECT_ID', 'CONSENT']
    mrnDF[consentFields].to_csv(filenames[0], sep="\t", header=True, index=False)

    # Sample Mapping
    ## What to use for the sample ID?  PhenotipsID + _1
    mrnDF['SAMPLE_ID'] = mrnDF['Phenotips_ID'] + '_1'
    ##SUE!! Remove Phenotips ID
    mapFields = ['Phenotips_ID', 'SUBJECT_ID', 'SAMPLE_ID']
    mrnDF[mapFields].to_csv(filenames[1], sep="\t", header=True, index=False)
    
    # Pedigree 
    mrnDF.rename(columns={'Phenotips_Family_ID': 'FAMILY_ID', 
                            'Mother_Phenotips_ID': 'MOTHER', 
                            'Father_Phenotips_ID': 'FATHER', 
                            'Gender' : 'SEX'}, inplace=True)
    pedFields = ['FAMILY_ID', 'SUBJECT_ID', 'MOTHER', 'FATHER', 'SEX']

    # fill in missing values for Mother and Father
    mrnDF.MOTHER.replace('^00$', '0', regex=True, inplace=True)
    mrnDF.MOTHER.replace(np.NaN, '0', inplace=True)
    mrnDF.MOTHER.replace("", '0', inplace=True)

    mrnDF.FATHER.replace(np.NaN, '0', inplace=True)
    mrnDF.FATHER.replace("", '0', inplace=True)

    # Convert gender to 1=Male 2=Female
    mrnDF.SEX.replace("M", '1', inplace=True)
    mrnDF.SEX.replace("F", '2', inplace=True)
    mrnDF.SEX.replace("", 'UNK', inplace=True)

    # Export pedigree file
    mrnDF[pedFields].to_csv(filenames[2], sep="\t", header=True, index=False)


    # Phenotypes
    ## gender affected, primary 0, no valid Age information in BSI?
    mrnDF.rename(columns={'Affected' : 'AFFECTION_STATUS', 
                            'Race' : 'RACE', 
                            'Ethnicity' : 'ETHNICITY'}, inplace=True)
    phenoFields = ['Phenotips_ID', 'AFFECTION_STATUS', 'SEX', 'RACE', 'ETHNICITY']

    # Replace to dictionary values
    mrnDF.AFFECTION_STATUS.replace("Y", '2', inplace=True)
    mrnDF.AFFECTION_STATUS.replace("N", '1', inplace=True)
    mrnDF.AFFECTION_STATUS.replace("U", 'UNK', inplace=True)

    ##SUE!! Are these MeSH
    mrnDF.RACE.replace("Am Indian/Alaska Nat", 'American Indian, North', inplace=True)
    mrnDF.RACE.replace("Asian", 'Asian', inplace=True)
    mrnDF.RACE.replace("Black/African Amer", 'African American', inplace=True)
    mrnDF.RACE.replace("Hawaiian/Pac. Island", 'Native Hawaiians', inplace=True)
    mrnDF.RACE.replace("Multiple Race", 'Multiple', inplace=True)
    mrnDF.RACE.replace("White", 'Caucasian', inplace=True)
    mrnDF.RACE.replace("Unknown", 'Unknown', inplace=True)

    ## SUE!! are these MeSH
    #mrnDF.ETHNICITY.replace("Not Hispanic", '1', inplace=True)
    #mrnDF.ETHNICITY.replace("Hispanic", '2', inplace=True)
    #mrnDF.ETHNICITY.replace("Unknown", '99', inplace=True)

    mrnDF[phenoFields].to_csv(filenames[3], sep="\t", header=True, index=False)

    # Sample Attributes
    ## gender affected, primary 0
    mrnDF['ANALYTE_TYPE'] = 'DNA'
    mrnDF['IS_TUMOR'] = 'no'
    mrnDF['Batch_Number'] = pd.to_numeric(mrnDF['Batch_Received'].str.replace("BATCH",""))

    # Batches 1-19?? are CIDR=1, Batches 20??-- are BCM=2
    mrnDF['SEQUENCING_CENTER'] = 'CIDR'
    mrnDF.loc[mrnDF.Batch_Number > 19, 'SEQUENCING_CENTER'] = 'BCM'

    mrnDF.rename(columns={'Tissue' : 'TISSUE'}, inplace=True)
    ##SUE!! using the coded values from BSI, MeSH would be better, but not all are clear values
    mrnDF.TISSUE.replace('9', 'Saliva', inplace=True)
    mrnDF.TISSUE.replace('15', 'Blood', inplace=True)
    mrnDF.TISSUE.replace('18', 'EBV Transformed B cell line', inplace=True)
    mrnDF.TISSUE.replace('21', 'PBMC', inplace=True)
    mrnDF.TISSUE.replace('22', 'Fibroblast', inplace=True)

    mrnDF['BODY_SITE'] = 'Unknown'
    mrnDF.loc[mrnDF.TISSUE == 'Saliva', 'BODY_SITE'] = 'Inner Oral Cavity'
    mrnDF.loc[mrnDF.TISSUE == 'Blood', 'BODY_SITE'] = 'Peripheral Blood'
    mrnDF.loc[mrnDF.TISSUE == 'EBV Transformed B cell line', 'BODY_SITE'] = 'Peripheral Blood'
    mrnDF.loc[mrnDF.TISSUE == 'PBMC', 'BODY_SITE'] = 'Peripheral Blood'
    mrnDF.loc[mrnDF.TISSUE == 'Fibroblast', 'BODY_SITE'] = 'Connective Tissue'

    # Export sample attribute file
    sampleFields = ['SAMPLE_ID', 'BODY_SITE', 'ANALYTE_TYPE', 'IS_TUMOR', 'TISSUE', 'SEQUENCING_CENTER']
    mrnDF[sampleFields].to_csv(filenames[4], sep="\t", header=True, index=False)

    # Done
    return()

####################################
# 
# Main 
#
####################################

def main():

    #
    # Usage statement
    #
    parseStr = 'Converts MRNs to Phenotips IDs\n\n\
    Usage:\n\
        hla_mrn2ptips.py -i infile \n\n\
    Example:\n\
        hla_mrn2ptips.py -i Lab78706121.csv \n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', required=False, nargs='?', type=argparse.FileType('r'), default=None, help='Input text file containing a list of data with MRNs to be converted to Phenotips IDs')
    parser.add_argument('-o', '--outfile', required=False, action='store', type=str, default="out.csv", help='Output file')
    parser.add_argument('-d', '--delim', required=False, type=str, default=",", help='column delimiter of infile, default=","')
    parser.add_argument('-c', '--column', required=False, type=int, default="1", help='column containing MRNs, default="1"')

    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    delim = args.delim
    colno = args.column

    #####################################
    #
    # Set up the log file
    #
    #####################################
    #thedate = str(datetime.datetime.now()).split()[0]
    #thedate = re.sub("-","",thedate)
    #global log
    #log = open('hla_mrn2ptips' + '.log', 'a')
    #log.write('\n' + str(datetime.datetime.now()) + '\n')
    #log.write(' '.join(sys.argv) + '\n')
    #log.write('hla_mrn2ptips.py version ' + __version__ + '\n\n')
    #log.flush()

    #####################################
    #
    # Read in the file and get the MRNs
    #
    #####################################
    test_file(infile)
    inDF = pd.read_csv(infile, header=0, encoding='utf-8')
    mrns = inDF['MRN'].astype(str)

    mrns = ["-".join([x[0:2], x[2:4], x[4:6], x[-1:]]) for x in mrns]
    inDF['MRN'] = mrns
    mrns = list(set(mrns))
    mrns.sort()

    ########################################
    #
    # Collection information from BSI, create the output files
    #
    ########################################
    # Set up the variables, bsi info
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "

    send_update("\nPulling sample data from BSI...")
    ## Establish a BSI Connection with user's credentials
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    ##
    ## Query BSI for all data through the current Batch
    ##
    # Pull data for subjects through the current batch
    fields = ['MRN', 'Phenotips ID', 'Batch Received']

    ##SUE!! Query fails with more than 615 MRNs, What if > 1000
    mrnDF = bsi_query(curl_get, url_reports, session, fields, mrns[:500], 'MRN')
    if len(mrns) > 500:
        mrnDF = pd.concat([mrnDF, bsi_query(curl_get, url_reports, session, fields, mrns[500:], 'MRN')])

    # Merge data with Phenotips IDs
    outDF = pd.merge(inDF, mrnDF, on='MRN', how='inner')
    outDF.drop(['MRN'], axis=1, inplace=True)

    # Move Phenotips ID to the beginning
    cols = list(outDF.columns.values)
    print(cols)
    cols = cols[-2:] + cols[:-2]
    outDF = outDF[cols]

    outDF.to_csv(outfile, header=True, index=False)


if __name__ == '__main__':
    main()

