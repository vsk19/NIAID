#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14, 2019
Updated on

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

csi_send_pedigree.py
    Generate a pedigree file to accompany CSI batches sent to CIDR for sequencing

v 1.0 - initial code version.

"""
__author__ = 'Susan Huse'
__date__ = 'February 14, 2019'
__version__ = '1.0'
__copyright__ = 'No copyright protection, can be used freely'

import sys
import os
import re
import datetime
import pandas as pd
import argparse
import numpy as np
#import docx2txt
#from docx import Document
from argparse import RawTextHelpFormatter
from ncbr_bsi import read_conf, send_curl, get_bsi_session, get_bsi_name, bsi_query, return_bsi_info
from ncbr_huse import send_update, err_out, pause_for_input, test_dir, test_file

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


# Read CIDR Manifest spreadsheet to get the Order IDs for the pedigree file
def pull_batch_ready():
    df = bsi_query(curl_get, url_reports, session, qfields, qcond, 'Batch Ready')
    return(df)
    
def get_all_family(familyIDs):
    df = bsi_query(curl_get, url_reports, session, qfields, familyIDs, 'Phenotips Family ID')
    return(df)

# Clean a few things up in the dataframe returned
def clean_bsi_data(df):

    # fill in missing values for Mother and Father
    df.Mother_Phenotips_ID.replace(np.NaN, '0', inplace=True)
    df.Mother_Phenotips_ID.replace("", '0', inplace=True)
    df.Mother_Phenotips_ID.replace('00', '0', inplace=True)
    df.Father_Phenotips_ID.replace(np.NaN, '0', inplace=True)
    df.Father_Phenotips_ID.replace("", '0', inplace=True)
    
    # Move gender word to letter
    df['Gender'].replace({'Female':'F', 'Male':'M'} ,inplace=True)
    df['Gender'] = np.where((df['Gender']=='F') | (df["Gender"]=='M'), df['Gender'], '0')

    
    # Remove family members that have already been sent
    df = df.drop(df[(df.Batch_Sent == "") & (~df.Batch_Ready.isin(qcond))].index)
    
    df = df.drop_duplicates()
    #    print("Order IDs:\n{}".format(df.head()))
    return(df)
    
def create_ped(df):
    #SUE!!
    #1) strip off FAM00* from Family ID -> Family
    #2) Individual if (proband 3, father 1, mother 2, sibling 4, etc)
    #3) Gender -> Sex
    #4) DNA_Class -> 1 for all valid samples, 0 otherwise
    #5) IRB -> NIAID
    #6) CRIS Order # -> Subject_ID
    #7) Race -> Population
    #8) Organism -> Human
    #9) Investigator_Column_1 -> relationship Mother, Father, Brother, Proband self,
    #10) Invesigator_Column_2 -> leave blank
    
    df = clean_bsi_data(df)
    
    # Remove prefix from Family
    df['Family'] = df['Phenotips_Family_ID'].str.replace('FAM0*', '', regex=True)
    
    # Add Columns
    #df['DNA_Class'] = np.where(df['Batch_Ready']=='CSI QC PASS', '1', '0')
    df['DNA_Class'] = 0
    df['IRB'] = "NIAID"
    df['Organism'] = 'Human'
    df['Investigator_Column_1'] = ''
    df['Investigator_Column_2'] = ''
    df['Father'] = 0
    df['Mother'] = 0
    
    # Set Proband
    # Proband Yes, Affected = Yes
    df['Individual'] = np.where((df['Proband']=='Yes') & (df['Affected'] == 'Y'), 3, 4)
    # Proband = Yes, and Relationship  = Proband
    df['Individual'] = np.where((df['Proband']=='Yes') & 
      (df["Relationship"].str.contains("Proband", case=False)), 3, df['Individual'])
    
    # Proband but Affected=No
    df['Individual'] = np.where((df['Proband']=='Yes') & (df['Affected'] != 'Y'), -3, df['Individual'])
    # Proband = No, but Relationship  = Proband
    df['Individual'] = np.where((df['Proband']!='Yes') & 
      (df["Relationship"].str.contains("Proband", case=False)), -3, df['Individual'])
    # Proband = Yes, but Relationship is not Proband
    df['Individual'] = np.where((df['Proband']=='Yes') & 
      ( (~df["Relationship"].str.contains("Proband", case=False)) & (df["Relationship"] != "")), -3, df['Individual'])

    # Set parents
    for m in df.Mother_Phenotips_ID:
        if m != '0':
            df['Individual'] = np.where((df['Phenotips_ID']== m) & (df['Relationship'].str.contains("Mother", case=False)), 
              2, df['Individual'])
            df['Individual'] = np.where((df['Phenotips_ID']== m) & (~df['Relationship'].str.contains("Mother", case=False)), 
              -2, df['Individual'])
            
            df['Mother'] = np.where(df['Mother_Phenotips_ID']== m, 2, df['Mother'])

    for f in df.Father_Phenotips_ID:
        if f != '0':
            df['Individual'] = np.where((df['Phenotips_ID']== f) & (df['Relationship'].str.contains("Father", case=False)), 
              1, df['Individual'])
            df['Individual'] = np.where((df['Phenotips_ID']== f) & (~df['Relationship'].str.contains("Father", case=False)), 
              -1, df['Individual'])
#            df['Individual'] = np.where(df['Phenotips_ID']== f, 1, df['Individual'])
            df['Father'] = np.where(df['Father_Phenotips_ID']== f, 1, df['Father'])
            
    # Check that there are no duplicate Individuals: multiple 4s
    singleton_families = []
    for f in df['Phenotips_Family_ID'].unique():
        x = df[(df['Phenotips_Family_ID'] == f) & (df['Individual'] > 3)]['Individual']
        if x.size > 1:
            x = 4
            for idx in df.index[(df['Phenotips_Family_ID'] == f) & (df['Individual'] > 3)].tolist():
                df.loc[idx, ['Individual']] = x
                x += 1
                
        # If the only family member at this point, call it a singleton
        if df[(df['Phenotips_Family_ID'] == f)].shape[0] == 1:
            singleton_families.append(f)

    # Check for missing parental Phenotips ID
    df['Individual'] = np.where((df['Relationship'].str.contains("Mother", case=False)) & (df['Individual'] != 2), -2, df['Individual'])
    df['Individual'] = np.where((df['Relationship'].str.contains("Father", case=False)) & (df['Individual'] != 1), -1, df['Individual'])
            
    df['Investigator_Column_1'] = np.where(df['Phenotips_Family_ID'].isin(singleton_families), "Singleton", df['Relationship'])

    # Final Clean up
    df.rename(index=str, inplace=True,
              columns={'Gender': 'Sex', 
                       'CRIS_Order#' : 'Subject_ID'})
    
#    print(df.to_string())    
    return(df.sort_values(['Phenotips_Family_ID']))
    
    
    
####################################
# 
# Main 
#
####################################

def main():

    #
    # Usage statement
    #
    parseStr = 'Creates a pedigree file to accompany sample batch sent for sequencing\n\
    Usage:\n\
        csi_send_pedigree.py -p pedfile -a allfile -f outformat\n\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-p', '--pedfile', required=False, action='store', type=str, 
                        default="Pedigree_for_CIDR.csv", help='Output file for pedigree information')
    parser.add_argument('-a', '--allfile', required=False, action='store', type=str, 
                        default="Pedigree_for_CIDR_all.csv", help='Output file with additional data supporting teh pedigree file')

    args = parser.parse_args()
#    infile = args.infile
#    batch = args.batch
    allfile = args.allfile
    pedfile = args.pedfile
    #outfile_path = "/Volumes/NIAID_Centralized_Sequencing/NCBR/to_CIDR"
    #outfile_ped_all = os.path.join(outfile_path, allfile)
    #outfile_ped = os.path.join(outfile_path, pedfile)

    ########################################
    #
    # Set up the variables, bsi info
    #
    ########################################
    global cnf, url_session, url_reports, curl_get, session, qfields, qcond
    
    cnf, url_session, url_reports, curl_get = return_bsi_info()
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    ########################################
    #
    # Collection information from BSI
    #
    ########################################
    qfields = ['Batch Ready', 'Phenotips Family ID', 'Phenotips ID', 
               'Father Phenotips ID', 'Mother Phenotips ID', 'Family Complete Status',
               'Proband', 'Affected', 'Gender',
               'CRIS Order #', 'CRIS Order Status', 'Patient Name',
               'Race', 'Relationship',
               'Batch Received', 'Batch Sent']
#    qcond = ['CSI QC PASS']
    qcond = ['Yes']
    
    df = pull_batch_ready()
    
    # Get full family membership
    familyIDs = df['Phenotips_Family_ID'].unique()
    familyDF = get_all_family(familyIDs)
   
    
    ########################################
    #
    # Collection information from BSI
    #
    ########################################
    pedDF = create_ped(familyDF)
    #print("Ped size: {}".format(pedDF.size()))
    #print("Returned Ped DF:\n{}".format(pedDF.head().to_string()))
    pedDF.to_csv(allfile, index=False)


    ######################################
    #   
    # Subset to Ped file columns only and write to the output file
    #
    ######################################
    final_fields = ['Family', 'Individual', 'Father', 'Mother', 'Sex', 
                    'DNA_Class', 'IRB', 'Subject_ID', 'Race',
                    'Organism', 'Investigator_Column_1', 'Investigator_Column_2']

    pedDF = pedDF[final_fields].to_csv(pedfile, index=False)

    ######################################
    #   
    # Close out and clean up
    #
    ######################################
#    send_update("\n{} successfully written by csi_send_pedigree.py".format(outfile))
    send_update(str(datetime.datetime.now()) + '\n')
    #log.close()

if __name__ == '__main__':
    main()

