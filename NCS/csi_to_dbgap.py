#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:16:58 2019
Updated on
Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical
csi_to_dbgap.py
    Reads the sample mapping and VCF information and 
    preps the metadata for uploading to GRIS
v 1.0 - initial code version.
    
"""


__author__ = 'Susan Huse'
__date__ = 'June 24, 2019'
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
from ncbr_huse import send_update, err_out
#import json
from ncbr_bsi import read_conf, get_bsi_session, bsi_query

######################################
#   
# Output the dbGaP Dataset (DS) files
#
######################################

def filter_results(df, batch):
    # only samples up through the specified batch number
    df = df[df['Batch_Received'] <= batch]

    # only data that are
    df = df[df['Batch_Received'].str.contains("BATCH")]

    # get rid of canceled orders
    df = df.loc[~df['CRIS_Order_Status'].str.contains('Cancel')]
    df = df[df['CRIS_Order_Status'] != 'Auto Complete']
    df.drop(columns = ['CRIS_Order_Status'], inplace=True)

    # require CIDR_Exome_ID SUE: will need to update this field to BCM
    df = df[df['CIDR_Exome_ID'] != '']
    
    # SUE check for unique CIDR Exome IDs on all 

    return(df)
    
def add_parent_ids(df):
    # pull mother and father IDs and look up their CIDR Exome IDs
#    print(df['Mother_Phenotips_ID'].head())
    mothers = list(set(df['Mother_Phenotips_ID']))
#    print(mothers[:10])
#    print(df['Father_Phenotips_ID'].head())
    fathers = list(set(df['Father_Phenotips_ID']))
#    print(fathers[:10])

    parent_fields = ['Phenotips ID', 'CIDR Exome ID', 'CRIS Order Status']
    mDF = bsi_query(curl_get, url_reports, session, parent_fields, mothers, 'Phenotips ID')
    mDF = mDF[mDF['CIDR_Exome_ID'] != '']
    mDF = mDF.loc[~mDF['CRIS_Order_Status'].str.contains('Cancel')]
    mDF = mDF[mDF['CRIS_Order_Status'] != 'Auto Complete']
    mDF.rename(columns={'Phenotips_ID': 'Mother_Phenotips_ID', 'CIDR_Exome_ID': 'MOTHER'}, inplace=True)
    mDF.drop(columns=['CRIS_Order_Status'], inplace=True)
    
    fDF = bsi_query(curl_get, url_reports, session, parent_fields, fathers, 'Phenotips ID')
    fDF = fDF[fDF['CIDR_Exome_ID'] != '']
    fDF = fDF.loc[~fDF['CRIS_Order_Status'].str.contains('Cancel')]
    fDF = fDF[fDF['CRIS_Order_Status'] != 'Auto Complete']
    fDF.drop(columns=['CRIS_Order_Status'], inplace=True)

    fDF.rename(columns={'Phenotips_ID': 'Father_Phenotips_ID', 'CIDR_Exome_ID': 'FATHER'}, inplace=True)
    
#    print(df.head())
#    print(mDF.head())
#    print(fDF.head())

    df = pd.merge(df, mDF, on='Mother_Phenotips_ID', how='left')
    df = pd.merge(df, fDF, on='Father_Phenotips_ID', how='left')
    
    # fill in missing values for Mother and Father
    df.MOTHER.replace('^00$', '0', regex=True, inplace=True)
    df.MOTHER.replace("0", '0', inplace=True)
    df.MOTHER.replace(np.NaN, '0', inplace=True)
    df.MOTHER.replace("", '0', inplace=True)

    df.FATHER.replace(np.NaN, '0', inplace=True)
    df.FATHER.replace("0", '0', inplace=True)
    df.FATHER.replace("", '0', inplace=True)

#    print(df.head())

    return(df)
    
def add_family_ids(df):
    # pull family IDs, get the proband, and use that CIDR Exome ID
#    df = dbgapDF
    probandDF = df[df['Proband'] == "Yes"][['Phenotips_Family_ID', 'Phenotips_ID', 'CIDR_Exome_ID']]
    probandDF['FAMILY_ID'] = "FAM-" + probandDF['CIDR_Exome_ID']
#    print(probandDF.head())

    df = pd.merge(df, probandDF, on='Phenotips_Family_ID', how='left')
    df.drop(columns=['Phenotips_ID_y', 'CIDR_Exome_ID_y'], inplace=True)
    df.rename(columns={'Phenotips_ID_x': 'Phenotips_ID', 'CIDR_Exome_ID_x': 'CIDR_Exome_ID'}, inplace=True)
    
    # There are a few that don't have family ID, remove them from the dataframe,
    # This is mostly due to proband not having been returned yet    
    lostDF = df[df.FAMILY_ID.isnull()]
    # print(lostDF)
    
    df = df[~df.FAMILY_ID.isnull()]
#    print(df.head())

    return(df)

######################################
#   
# Output the dbGaP Dataset (DS) files
#
######################################

def write_files(df, filenames):
    #filenames 0: file_consent, 1:file_mapping, 2:file_ped, 3:file_pheno, 4:file_samples):

    #
    # Consent - ID and Consent = 1
    #
    ## SUE: we need to be sure that if we have more than one sample we use the consistent Subject ID
    df['SUBJECT_ID'] = df['CIDR_Exome_ID']
    df['CONSENT'] = 1
    ###dbgapDF['SUBJECT_SOURCE'] = 'NIAID CSI'
    
    ###dbgapDF['SOURCE_SUBJECT_ID'] = dbgapDF['Phenotips_ID']
    ###consentFields = ['SUBJECT_ID', 'CONSENT', 'SUBJECT_SOURCE', 'SOURCE_SUBJECT_ID']
    consentFields = ['SUBJECT_ID', 'CONSENT']
    df[consentFields].to_csv(filenames[0], sep="\t", header=True, index=False)

    #
    # Sample Mapping
    #
    ## What to use for the sample ID?  PhenotipsID + _1
    ## SUE: we need to be sure that if we have more than one sample we use the consistent
    ## Subject ID and increment the sample ID
    df['SAMPLE_ID'] = df['SUBJECT_ID'] + '_1'
    mapFields = ['SUBJECT_ID', 'SAMPLE_ID']
    df[mapFields].to_csv(filenames[1], sep="\t", header=True, index=False)
    
    #
    # Pedigree 
    #
    #http://www.gwaspi.org/?page_id=145
    #The fields in a PED file are
    #Family ID
    #Sample ID
    #Paternal ID
    #Maternal ID
    #Sex (1=male; 2=female; other=unknown)
    #Affection (0=unknown; 1=unaffected; 2=affected)
    #Genotypes (space or tab separated, 2 for each marker. 0=missing)

    df.rename(columns={'Gender' : 'SEX'}, inplace=True)
    pedFields = ['FAMILY_ID', 'SUBJECT_ID', 'MOTHER', 'FATHER', 'SEX']


    # Convert gender to 1=Male 2=Female
    df.SEX.replace("M", '1', inplace=True)
    df.SEX.replace("F", '2', inplace=True)
    df.SEX.replace("", 'UNK', inplace=True)


    #
    # Export pedigree file
    #
    df[pedFields].to_csv(filenames[2], sep="\t", header=True, index=False)


    #
    # Phenotypes
    #
    ## gender affected, primary 0, no valid Age information in BSI?
    df.rename(columns={'Affected' : 'AFFECTION_STATUS', 
                            'Race' : 'RACE', 
                            'Ethnicity' : 'ETHNICITY'}, inplace=True)
    phenoFields = ['SUBJECT_ID', 'AFFECTION_STATUS', 'SEX', 'RACE', 'ETHNICITY']

    # Replace to dictionary values
    df.AFFECTION_STATUS.replace("Y", '2', inplace=True)
    df.AFFECTION_STATUS.replace("N", '1', inplace=True)
    df.AFFECTION_STATUS.replace("U", 'UNK', inplace=True)

    ##SUE!! Are these MeSH
    df.RACE.replace("Am Indian/Alaska Nat", 'American Indian, North', inplace=True)
    df.RACE.replace("Asian", 'Asian', inplace=True)
    df.RACE.replace("Black/African Amer", 'African American', inplace=True)
    df.RACE.replace("Hawaiian/Pac. Island", 'Native Hawaiians', inplace=True)
    df.RACE.replace("Multiple Race", 'Multiple', inplace=True)
    df.RACE.replace("White", 'Caucasian', inplace=True)
    df.RACE.replace("Unknown", 'Unknown', inplace=True)

    ## SUE!! are these MeSH
    #dbgapDF.ETHNICITY.replace("Not Hispanic", '1', inplace=True)
    #dbgapDF.ETHNICITY.replace("Hispanic", '2', inplace=True)
    #dbgapDF.ETHNICITY.replace("Unknown", '99', inplace=True)

    df[phenoFields].to_csv(filenames[3], sep="\t", header=True, index=False)

    #
    # Sample Attributes
    #
    ## gender affected, primary 0
    df['ANALYTE_TYPE'] = 'DNA'
    df['IS_TUMOR'] = 'no'
    df['Batch_Number'] = pd.to_numeric(df['Batch_Received'].str.replace("BATCH",""))

    # Batches 1-19?? are CIDR=1, Batches 20??-- are BCM=2
    df['SEQUENCING_CENTER'] = 'CIDR'
    df.loc[df.Batch_Number > 19, 'SEQUENCING_CENTER'] = 'BCM'

    df.rename(columns={'Tissue' : 'BODY_SITE'}, inplace=True)
    ##SUE!! using the coded values from BSI, MeSH would be better, but not all are clear values
    df.BODY_SITE.replace('9', 'Saliva', inplace=True)
    df.BODY_SITE.replace('15', 'Blood', inplace=True)
    df.BODY_SITE.replace('18', 'EBV Transformed B cell line', inplace=True)
    df.BODY_SITE.replace('21', 'PBMC', inplace=True)
    df.BODY_SITE.replace('22', 'Neutrophils', inplace=True)
    df.BODY_SITE.replace('23', 'Fibroblast', inplace=True)

    ##SUE: get rid of body site, just use tissue
    #df['BODY_SITE'] = 'Unknown'
#    df.loc[df.TISSUE == 'Saliva', 'BODY_SITE'] = 'Inner Oral Cavity'
#    df.loc[df.TISSUE == 'Blood', 'BODY_SITE'] = 'Peripheral Blood'
#    df.loc[df.TISSUE == 'EBV Transformed B cell line', 'BODY_SITE'] = 'Peripheral Blood'
#    df.loc[df.TISSUE == 'PBMC', 'BODY_SITE'] = 'Peripheral Blood'
#    df.loc[df.TISSUE == 'Fibroblast', 'BODY_SITE'] = 'Connective Tissue'

    # Export sample attribute file
    sampleFields = ['SAMPLE_ID', 'ANALYTE_TYPE', 'IS_TUMOR', 'BODY_SITE', 'SEQUENCING_CENTER']
    df[sampleFields].to_csv(filenames[4], sep="\t", header=True, index=False)

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
    parseStr = 'Creates dbGaP Dataset files for all available data through the specified Batch Received #.\n\n\
    Usage:\n\
        csi_to_dbgap.py -b batch \n\n\
    Example:\n\
        csi_to_dbgap.py -b 7\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-b', '--batch', required=True, type=int, help='Batch number (integer)')
    #parser.add_argument('-p', '--pedfile_only', required=False, action='store_true', default=False, help='Output ped file only, do not update the tracking, masterkey or other files')
    #parser.add_argument('-d', '--dir_info', required=False, type=str, default="./rawdata/Sample_Info", help='Directory containing input CIDR csv files ("rawdata/Sample_Info/")')

    args = parser.parse_args()
    batch = args.batch
    pd.set_option('display.max_columns', None)
    
    #####################################
    #
    # Set up the variables and the log file
    #
    #####################################

    # Set up the log file
    thedate = str(datetime.datetime.now()).split()[0]
    thedate = re.sub("-","",thedate)
    global log 
    log = open('csi_to_dbgap' + '.log', 'a')
    log.write('\n' + str(datetime.datetime.now()) + '\n')
    log.write(' '.join(sys.argv) + '\n')
    log.write('csi_to_dbgap.py version ' + __version__ + '\n\n')
    log.flush()

    #####################################
    #
    # Load the Config File, and check it is complete
    #
    #####################################
    ## Set Batch name (BSI value) and Batch label (output text value)
    if batch < 10: 
       batch_name = "BATCH0" + str(batch)
    else:
       batch_name = "BATCH" + str(batch)
    #batch_label = "Batch " + str(batch)

    # Set up the variables, bsi info
    global url_session, url_reports, curl_get, cnf, user, pw, session
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    ########################################
    #
    # Collection information from BSI, create the output files
    #
    ########################################
    send_update("\nPulling sample data from BSI...")
    ## Establish a BSI Connection with user's credentials

    ##
    ## Query BSI for all data through the current Batch
    ##
    # Pull data for subjects through the current batch
    fields = ['Phenotips ID', 
              'Phenotips Family ID', 
              'Father PhenotipsId', 
              'Mother PhenotipsId', 
              'CIDR Exome ID',
              'Gender', 
              'Proband',
              'Affected Status', 
              'Batch Received', 
              'Race',
              'Ethnicity',
              'Age', 
              'Tissue Origin',
              'CRIS Order Status', 
              'Order Date']

    dbgapDF = bsi_query(curl_get, url_reports, session, fields, ['BATCH*'], 'Batch Received', False)
#    print(dbgapDF.shape)

    dbgapDF = filter_results(dbgapDF, batch_name)
#    print(dbgapDF.shape)
    
    dbgapDF = add_parent_ids(dbgapDF)
#    print(dbgapDF.head())
    
    dbgapDF = add_family_ids(dbgapDF)
#    print(dbgapDF.head())

    ######################################
    #   
    # Creating output files (SampleTracking, Masterkey, Pedigree, Batch Info
    #
    ######################################
    send_update("\nWriting output files...", log)

    # Print MasterKey, BatchInfo, Pedigree files
    ##SUE: add back in the 2a. and numeric prefix to file names
    outfiles = ['2a_SubjectConsent_DS.txt', 
                '3a_SSM_DS.txt', 
                '4a_Pedigree_DS.txt', 
                '5a_SubjectPhenotypes_DS.txt', 
                '6a_SampleAttributes_DS.txt']

    write_files(dbgapDF, outfiles)

    ######################################
    #   
    # Close out and clean up
    #
    ######################################
    send_update("\ncsi_to_dbgap.py successfully completed", log)
    send_update(str(datetime.datetime.now()) + '\n', log)
    log.close()

if __name__ == '__main__':
    main()