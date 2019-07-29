#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9, 2019
Updated on

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

CSI_neg_rpts.py
    Generates csv file with info on patients for which a negative molecular diagnostic report is to be written

v 1.0 - initial code version.

"""
__author__ = 'Susan Huse'
__date__ = 'January 9, 2019'
__version__ = '1.0'
__copyright__ = 'No copyright protection, can be used freely'

import sys
import os
import re
import datetime
import pandas as pd
import argparse
import datetime
import time
from argparse import RawTextHelpFormatter
from ncbr_bsi import read_conf, send_curl, get_bsi_session, get_bsi_name, bsi_query, return_bsi_info
from ncbr_huse import send_update, test_dir, test_file

####################################
# 
# Functions for processing data
#
####################################
# Test file exists or error out
def write_out(df, fn, fmt):
    if fmt == "csv":
        df.to_csv(fn, index=False) #, quoting=1)
    else:
        df.to_excel(fn, sheet_name='negatives', index=False)

# Read in the MRNs - clean and create a df of MRN:Test
def read_mrn_test(infile):
    lines = [line.strip() for line in infile]
    mrns = []
    tests = []
    testname="Uknown"
    comments = []
    comment = ""

    # For each line sort as MRN, Batch# or test name
    rMRN = re.compile('\d{2}-\d{2}-\d{2}-\d')
    rBatch = re.compile('Batch')

    for line in lines:
        # skip blank lines
        if line == "":
            next

        if line.isspace():
            next

        # skip batch numbers
        elif rBatch.match(line) is not None:
            next

        # remove trailing colons
        if line.endswith(':'):
            line = line[:-1]

        # if it is an MRN, pull it, and set to the last test name read
        # if MRN put all trailing text into a Comments field
        if rMRN.match(line) is not None:
            mrns.append(rMRN.search(line).group())
            tests.append(testname)

            #print(", ".join(mrns))
            #print(", ".join(tests))
            comments.append(re.sub(mrns[-1], "", line).strip())

        # otherwise must be test name
        else:
            testname = line
    
    # Create a df of all the MRNs and test names
    #testDF = pd.DataFrame.from_items([('MRN', mrns), ('Test', tests)])
    testDF = pd.DataFrame({'MRN': mrns, 'Test':tests, 'TestComment':comments})
    #print(testDF)
    return(testDF)

def get_specimen_type(row):
    if row['TestName'] == "Centralized Exome, Blood, NIAID":
        return("Blood")
    if row['TestName'] == "Centralized Exome, Saliva, NIAID":
        return("Saliva")
    if row['TestName'] == "Centralized Exome, DNA, NIAID":
        return("noted in the SpecialInstructions")
    if row['TestName'] == "Centralized Exome Order, NIAID":
        return("noted in the SpecialInstructions")
    if row['TestName'] == "Extracted DNA":
        return("Unknown")

def make_date_type(x):
    if(type(x) == 'datetime.date' or type(x) != 'str'):
        x = datetime.date.strftime(x, "%Y-%m-%d %H:%M")
    x = datetime.datetime.strptime(x, '%Y-%m-%d %H:%M').date()
    return(x)

def get_collect_date(row):
    oldDate = datetime.datetime.strptime('10012018', '%m%d%Y').date()
    newDate = make_date_type(row['OrderDate'])
    #newDate = datetime.date.strftime(row['OrderDate'], "%Y-%m-%d %H:%M")
    #newDate = datetime.datetime.strptime(newDate, '%Y-%m-%d %H:%M').date()
    rCollectDate = re.compile('date of collection')

    if newDate < oldDate:
        if rCollectDate.match(row['SpecialInstructions'], re.IGNORECASE) is not None:
            x = re.sub('date of collection:?', '', row['SpecialInstructions'])
            print(row['SpecialInstructions'])
            print(x)
    if row['Specimen Type'] == "Blood":
        return(row['BloodCollectDate'])
    if row['Specimen Type'] == "Saliva":
        return(row['OriginalSpecimenCollectionDate'])

# Pull data from GRIS daily dumps
def read_CRIS_dump(gris_dir, mrns, cris_fields, status_values): 
    # Get the list of available dumps, and grab the last one
    dump_names = os.listdir(gris_dir)
    dump_names = [f for f in dump_names if re.search(".xls",f)]
    dump_names.sort()
    dump_file = os.path.join(gris_dir, dump_names[-1])

    # Read in the file
    df = pd.read_excel(dump_file, sheet_name="RawData")

    # Filter on CRIS fields, MRNs, and Order Status
    df = df[cris_fields]
    df = df.loc[df['MRN'].isin(mrns)]
    df = df[df['OrderStatus'].isin(status_values)]

    # remove newlines and commas that wreak havoc with csv files
    df['SpecialInstructions'] = df['SpecialInstructions'].str.replace('\n',' ')
    df['SpecialInstructions'] = df['SpecialInstructions'].str.replace(',',';')

    # Add SpecimenType field as: Blood, Saliva, or "Noted in SpecialInstructions"
    df['Specimen Type'] = df.apply(lambda row: get_specimen_type(row), axis=1)
    #df['DATE_OF_COLLECTION'] = df.apply(lambda row: get_collect_date(row), axis=1)

    return(df)

# Pull data from BSI
def query_BSI_data(mrns, bsi_fields, status_values):

    cnf, url_session, url_reports, curl_get = return_bsi_info()
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    # Return all records for the MRNs, with correct Order Status
    df = bsi_query(curl_get, url_reports, session, bsi_fields, mrns, 'MRN')
    df = df[df['CRIS_Order_Status'].isin(status_values)]

    return(df)


####################################
# 
# Main 
#
####################################

def main():

    #
    # Usage statement
    #
    parseStr = 'Creates a CSV file with required fields for GRIS to create a series of negatives reports\n\
    from a text file containing MRNs associated with specific negative tests.\n\
    Usage:\n\
        CSI_neg_rpts.py -i infile -o outfile -f outformat\n\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', required=True, nargs='?', type=argparse.FileType('r'), default=None, 
                        help='Input text file containing a list of MRNs\nfor which negative reports are to be generated')
    parser.add_argument('-o', '--outfile', required=False, action='store', type=str, 
                        default="Negative_Reports", help='Output file for CSI report')
    #parser.add_argument('-t', '--test', required=False, action='store', type=str, 
    #                    help='Name of genetic test applied that is to be included in the negative report')
    parser.add_argument('-d', '--gris_dir', required=False, action='store', type=str,
                        default="/Volumes/NIAID_Centralized_Sequencing/GRIS/CRIS_DATA_DUMP/EXO",
                        help='Directory containing GRIS "CentralizedExomeOrders_*.xls"\n\
                        (default="/Volumes/NIAID_Centralized_Sequencing/GRIS/CRIS_DATA_DUMP/EXO")')
    parser.add_argument('-f', '--out_format', required=False, type=str, default='csv',
                        help='Output file format: csv or xlsx, default=csv', choices=['csv','xlsx'])
    args = parser.parse_args()
    infile = args.infile
    gris_dir = args.gris_dir

    outfile = args.outfile
    out_format = args.out_format
    outfile = outfile + "." + out_format

    #####################################
    #
    # Set up the log file
    #
    #####################################
    ## Set up the log file
    #thedate = str(datetime.datetime.now()).split()[0]
    #thedate = re.sub("-","",thedate)
    #global log 

    #log = open('csi_reporter' + '.log', 'a')
    #log.write('\n' + str(datetime.datetime.now()) + '\n')
    #log.write(' '.join(sys.argv) + '\n')
    #log.write('csi_reporter.py version ' + __version__ + '\n\n')
    #log.flush()

    #
    # Test the input file and directory
    #
    if infile:
        test_file(infile)
    
    test_dir(gris_dir)

    # 
    # List all the fields and filter values
    #
    cris_fields = ['MRN', 
                   'ClientDisplayName',
                   'CRISOrderID', 
                   'SpecimenCollected', 
                   'BloodCollectDate', 
                   'OriginalSpecimenCollectionDate', 
                   'SpecialInstructions', 
                   'TestName', 
                   'OrderDate',
                   'OrderStatus']

    status_values = ['Specimen Collected', 'Auto Complete']

    bsi_fields = ['MRN',
              'CRIS Order Status',
              'GRIS Owner', 
              'Batch Received',
              'CIDR Exome ID']

    #
    # Read in the MRNs - clean and create a df of MRN:Test
    #
    testDF = read_mrn_test(infile)
    mrns = testDF['MRN'].tolist()
    print(testDF.head())

    #   
    # Pull data from GRIS daily dumps
    #   
    dumpDF = read_CRIS_dump(gris_dir, mrns, cris_fields, status_values)
    #print("CRIS Data Dump:\n{}".format(dumpDF.head()))

    #####################################
    #
    # MRNs need to be quoted to get the hyphens through the curl command
    # Also BSI -> Tools -> User Profile -> Preferences -> Reports -> Range Delimiter: "~" not "-"
    #
    #####################################
    mrns = ["'"+x+"'" for x in mrns]
    #print("MRNs:\n{}".format(mrns))

    ########################################
    #
    # Collection information from BSI
    #
    ########################################
    bsiDF = query_BSI_data(mrns, bsi_fields, status_values)
    #print("BSI:\n{}".format(bsiDF.head()))

    ######################################
    #   
    # Merge the two dataframe sources, add test and reorder
    #
    ######################################
    df = pd.merge(dumpDF, bsiDF, on="MRN")
    df = pd.merge(df, testDF, on="MRN")

    ######################################
    #   
    # Filter data, check for dupes and drops, fix date formats
    #
    ######################################

    # Check for duplicate MRNs
    dupeMRNs = df['MRN'][df['MRN'].duplicated()]
    if len(dupeMRNs.shape) > 0:
        print("\n\nDuplicate MRNs returned in query:\n{}".format(", ".join(dupeMRNs.tolist())))

    # remove time from dates so loads correctly into excel
    date_fields = ['SpecimenCollected'] #, 'BloodCollectDate', 'OriginalSpecimenCollectionDate']
    for c in date_fields:
        #dates = [str(x).split(" ")[0] for x in df[c]]
        #dates = ["" if x == "NaT" else x for x in dates]
        dates = [make_date_type(x) for x in df[c]]
        df[c] = dates
    
    ######################################
    #   
    # Reset column names and write to the output file
    #
    ######################################
    df['Comments'] = ""
    df = df[['ClientDisplayName',
             'MRN',
             'CRISOrderID', 
             'BloodCollectDate',
             'OriginalSpecimenCollectionDate',
             'SpecimenCollected', 
             #'DATE_OF_COLLECTION',
             'Test',
             'GRIS_Owner',
             'Batch_Received',
             'Specimen Type', 
             'SpecialInstructions', 
             'TestComment',
             'Comments']]

    new_fields = {'ClientDisplayName':'Patient Name',
                  'MRN':'MRN',
                  'CRISOrderID':'Order_ID', 
                  'BloodCollectDate':'BloodCollectDate',
                  'OriginalSpecimenCollectionDate':'OriginalSpecimenCollectionDate',
                  'SpecimenCollected':'DATE_RECEIVED',
                  #'DATE_OF_COLLECTION':'DATE_OF_COLLECTION',
                  'Test':'TEST',
                  'GRIS_Owner':'GRIS_Owner',
                  'Batch_Received':'Batch Received',
                  'Specimen Type':'Specimen Type', 
                  'SpecialInstructions':'SpecialInstructions',
                  'Test_Comments':'TestComment',
                  'Comments':'Comments'}

    df.rename(index=str, columns=new_fields)

    write_out(df, outfile, out_format)

    ######################################
    #   
    # Close out and clean up
    #
    ######################################
    print("\n{} successfully written by CSI_neg_rpts.py".format(outfile))
    #send_update("\n{} successfully written by CSI_neg_rpts.py".format(outfile))

if __name__ == '__main__':
    main()

