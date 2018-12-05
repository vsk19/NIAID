#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:10:09 2018

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

ncs_to_bsi.py
    takes lists of ids and fields and queries BSI

v 1.0 - initial code version.

"""
__author__ = 'Susan Huse'
__date__ = 'December 4, 2018'
__version__ = '1.0'
__copyright__ = 'No copyright protection, can be used freely'

import sys
import os
import re
import datetime
import pandas as pd
import numpy as np
from yaml import load
import argparse
from argparse import RawTextHelpFormatter
from ncbr_huse import send_update, err_out, pause_for_input

import requests
import requests.auth
from requests.auth import HTTPDigestAuth

####################################
# 
# Functions for processing data
#
####################################
# Test file exists or error out
def test_file(f):
    if not os.path.isfile(f):
        err_out("Error: unable to locate input file:  " + f, log)
    else:
        return(True)

# Create the prompt to be used by get_bsi_data
def bsi_prompt_text(bsifile, idname, reccnt):
    return('\n'.join(["", stars, stars, "Please log into BSI and create a new report that includes:",
    "Subject ID", "PhenotipsId", "CRIS Order #", "Phenotips Family ID", "Batch Sent\n",
    "Use the following {}s as search criteria ({} =@(like) *inserted IDs*).".format(idname, idname),
    'Run the report and save the file as: "{}" in the same directory as {}'.format(bsifile, config['config_fname']),
    '\nPlease enter "y" to print out the list of {} {}s to use, or "q" to quit.\n'.format(reccnt, idname)]))

# Get the Phenotips information for the subjects in this batch and previous batches
def get_bsi_data(f1, f2, mapidx, added, dropped):
    # Prompt and import the first set of BSI data - subjects released in this batch
    # Get the full pedigree for all of the samples, received, added, and dropped
    print("MapIDx type={}:\n{}".format(type(mapidx), mapidx))
    print("Dropped type={}:\n{}".format(type(dropped), dropped))
    ## SUE: unhashable type: list
    rec_sent = list(set([mapidx, dropped]))
    print("Received and Sent {}:\n{}",format(str(len(rec_sent),rec_sent)))

    if not testmode:
        pause_for_input(bsi_prompt_text(f1, "CRIS Order #", str(len(rec_sent))), 'y', 'q', log)
        print("\n" + "\n".join(rec_sent, ))
        pause_for_input('\nPlease enter "y" to continue when you have created {}, or "q" to quit.\n'.format(f1), 'y','q', log)

    bsi1 = import_bsi_phenotips(f1)
    print("BSI1 Columns:\n{}".format(bsi1.columns.tolist()))
    print("BSI1:\n{}".format(bsi1.head()))
    
    # Pull all family members that have been sent
    if not testmode:
        familyIDs = list(set(bsi1['Phenotips_Family_ID']))
        pause_for_input(bsi_prompt_text(f2, "Phenotips Family ID", str(len(familyIDs))), 'y', 'q', log)
        print("\n" + "\n".join(familyIDs))
        pause_for_input('\nPlease enter "y" to continue when you have created this file, {}, or "q" to quit.\n'.format(f2), 'y','q', log)

    bsi2 = import_bsi_phenotips(f2)
    print("BSI2 Columns:\n{}".format(bsi2.columns.tolist()))
    print("BSI2:\n{}".format(bsi2.head()))
    
    bsi = pd.concat([bsi1,bsi2], sort=True)
    bsi.drop_duplicates(inplace=True)
    bsi = bsi.set_index(['Subject_ID'])
    #print("BSI:\n{}".format(bsi.head()))
    return(bsi)

# Create the prompt to be used by get_bsi_exomes
def bsi_exomes_prompt_text(bsifile, idname, reccnt):
    return('\n'.join(["", stars, stars, "Please log into BSI and create a new report that includes:",
    "CRIS Order #", "CIDR Exome ID", "Batch Received\n",
    "Use the following {}s as search criteria ({} =@(like) *inserted IDs*).".format(idname, idname),
    'Run the report and save the file as: "{}" in the same directory as {}'.format(bsifile, config['config_fname']),
    '\nPlease enter "y" to print out the list of {} {}s to use, or "q" to quit.\n'.format(reccnt, idname)]))

# Get the CIDR_Exome_ID information for the subjects missing exome IDs
def get_bsi_exomes(f3, ids):
    # Prompt and import the set of missing data 
    if not testmode:
        pause_for_input(bsi_exomes_prompt_text(f3, "CRIS Order #", str(len(ids))), 'y', 'q', log)
        print("\n" + "\n".join(ids))
        pause_for_input('\nPlease enter "y" to continue when you have created {}, or "q" to quit.\n'.format(f3), 'y','q', log)

    bsi = import_bsi_exomes(f3)
    bsi = bsi.set_index(['Subject_ID'])
    bsi = bsi.rename(columns = { 'PhenotipsId':'Phenotips_ID', 
                                 'Phenotips Family ID': 'Phenotips_Family_ID',
                                 'Batch Sent':'Batch_Sent', 
                                 'CIDR Exome ID':'CIDR_Exome_ID',
                                 'Father PhenotipsID':'Father_Phenotips_ID',
                                 'Mother PhenotipsID':'Mother_Phenotips_ID'})
    return(bsi)

######################################
#   
# Functions for importing data
#
######################################

# Import the BSI phenotips information
def import_bsi_phenotips(f):
    test_file(f)

    # import and convert to dataframe, using field names
    phenotips = pd.read_csv(f, header=3)

    # Test that all necessary columns are there
    theCols = phenotips.columns
    missing = []

    for c in ['Subject ID', 'PhenotipsId', 'CRIS Order #', 'Phenotips Family ID', 'Batch Sent', 'Father PhenotipsId', 'Mother PhenotipsId', 'Gender', 'Proband']:
        if c not in theCols:
            missing.append(c)

    if len(missing) > 0:
        err_out("Phenotips CSV file exported from BSI is missing required field(s):\n{}".format("\n".join(missing)), log)
    
    # rename the columns
    phenotips = phenotips.rename(columns = {'Subject ID':'MRN', 
                                            'PhenotipsId':'Phenotips_ID', 
                                            'Phenotips Family ID': 'Phenotips_Family_ID',
                                            'CRIS Order #':'Subject_ID', 
                                            'CIDR Exome ID':'CIDR_Exome_ID',
                                            'Batch Sent':'Batch_Sent', 
                                            'Father PhenotipsId':'Father_Phenotips_ID',
                                            'Mother PhenotipsId':'Mother_Phenotips_ID'})
    
    # drop the row count column
    phenotips = phenotips.drop(columns=['OBS'])

    return(phenotips)
    
# Import the BSI CIDR Exome ID information
def import_bsi_exomes(f):
    test_file(f)

    # import and convert to dataframe, using field names
    exomes = pd.read_csv(f, header=2)

    # Test that all necessary columns are there
    theCols = exomes.columns
    missing = []

    for c in ['CRIS Order #', 'CIDR Exome ID', 'Batch Received']:
        if c not in theCols:
            missing.append(c)

    if len(missing) > 0:
        err_out("CIDR Exome ID CSV file (bsi3.csv) exported from BSI is missing required field(s):\n{}".format("\n".join(missing)), log)
    
    # rename the columns
    exomes = exomes.rename(columns = {'CRIS Order #':'Subject_ID', 
                                            'CIDR Exome ID':'CIDR_Exome_ID',
                                            'Batch Received':'Batch_Received'})
    
    # drop the row count column
    exomes = exomes.drop(columns=['OBS'])

    return(exomes)
    
####################################
# 
# Main 
#
####################################

def main():
    # Global variables
    global batch
    global testmode
    global testurl

    #
    # Usage statement
    #
    parseStr = 'Reads a list of fields and IDs and queries BSI,\n\
    Usage:\n\
        ncs_to_bsi.py -i ids -f fields \n\n\
    Example:\n\
        ncs_to_bsi.py -b 7 -c batch7.config\n'

    
    #parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    #parser.add_argument('-b', '--batch', required=True, type=str, help='Batch ID, e.g., BATCH09 or BATCH12')
    #parser.add_argument('-i', '--ids', required=True, type=str, help='List of IDs')
    #parser.add_argument('-q', '--queryfield', required=True, type=str, help='filename to search query')
    #parser.add_argument('-f', '--fields', required=False, action='store_true', default=False, help='list of fields to include in query')

    #args = parser.parse_args()
    #batch = args.batch
    #ids = args.ids
    #queryField = args.queryfield
    #fields = args.fields
    #testmode=args.testmode

    #####################################
    #
    # Set up the variables and the log file
    #
    #####################################

    ## Set up the log file
    #thedate = str(datetime.datetime.now()).split()[0]
    #thedate = re.sub("-","",thedate)
    #global log 
    #log = open('ncs_to_bsi' + '.log', 'a')
    #log.write('\n' + str(datetime.datetime.now()) + '\n')
    #log.write(' '.join(sys.argv) + '\n')
    #log.write('ncs_to_bsi.py version ' + __version__ + '\n\n')
    #log.flush()

    #####################################
    ##
    ## Import each of the reference files, do pedigree first to find the duplicate sample
    ##
    #####################################
    cnf = "/Users/husesm/.my.cnf"
    batch = "BATCH10"
    ids = [ '002FSSVSS', '002FWJVNQ', '002FWJVPB', '002FWSBQL', '002FWSBQN', '002FWWSNV', '002FWZYHM', '002FXGWGF', '002FXLBHB', '002FXLRVP', '002FXNMHK', '002FXNMHV', '002FXNMRH', '002FXQBPC', '002FXQBVQ', '002FXQBWD', '002FXQBYN', '002FXQCNK', '002FXQSQL', '002FXZGBY', '002FXZGHP']

    fields = ['CRIS Order #', 'Batch Sent', 'Batch Received', 'CRIS Order Status', 'Subject ID', 'PhenotipsId', 'Phenotips Family ID']
    queryField = 'CRIS Order #'
    #print ("Batch, QueryField:{}, {}".format(batch, queryField))
    
    testurl = "https://rest.bsisystems.com/doc/dist/?url=https://rest.bsisystems.com/api/rest/swagger.json#/common"

    curl = "curl -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: EBMS+P9M&#33;76irHGx2q1a9L$n6a*1Hmjl' 'https://rest.bsisystems.com/api/rest/EBMS/reports/list?display_fields=subject.subject_id&display_fields=sample.field_252&display_fields=sample.field_274&display_fields=subject_131.field_170&criteria=subject.study_id%3DNIAID%20Centralized%20Sequencing&limit=10000&type=1'"
    url = "https://rest.bsisystems.com/api/rest/EBMS/reports/list?display_fields=sample.field_252&display_fields=sample.field_274&display_fields=subject_131.field_170&criteria=subject.study_id%3DNIAID%20Centralized%20Sequencing&limit=10000&type=1"
    # subject id, 252= Phenotips ID, 274= Cris order #, 170=Phenotiups Family ID, 

    def get_bsi_session():
        print("Get Session ID")
        url = 'https://httpbin.org/digest-auth/auth/user/pass'
        requests.get(url, auth=HTTPDigestAuth('user', 'pass'))
    
    def read_conf(cnf):
        with open(cnf, 'r') as f:
            x = f.readlines()

        user = x[0].rstrip()
        pw = x[1].rstrip()
        return(user, pw)

    user, pw = read_conf(cnf)
    #print("User {} and pw {}".format(user, pw))
    #session = get_bsi_session()
    
    #url = 'https://www.googleapis.com/qpxExpress/v1/trips/search?key=mykeyhere'
    #payload = open("request.json")
    #headers = {'content-type': 'application/json', 'Accept-Charset': 'UTF-8'}
    r = requests.get(url, auth=HTTPDigestAuth(user, pw))
    
    print("Request status code is OK: {}".format(r.status_code == requests.codes.ok))

    ######################################
    #   
    # Close out and clean up
    #
    ######################################
    #send_update("\nncs_to_bsi.py successfully completed", log)
    #send_update(str(datetime.datetime.now()) + '\n', log)
    #log.close()

if __name__ == '__main__':
    main()



