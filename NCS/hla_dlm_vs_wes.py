#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 08:50:58 2019
Updated on

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

hla_mrn2ptips.py
    Compare the CSI pipeline HLA typing with the CC typing

v 1.0 - initial code version.
    
"""
__author__ = 'Susan Huse'
__date__ = 'June 28, 2019'
__version__ = '1.0'
__copyright__ = 'No copyright protection, can be used freely'

#import sys
import os
#import re
#import datetime
import pandas as pd
#import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from ncbr_huse import test_file
#import fnmatch

######################################
#   
# Functions
#
######################################

# Read in the CC test values 
def read_CC(cc_file):
#    test_file(cc_file)
    df = pd.read_csv(cc_file, header=0, encoding='utf-8')
    df = df[ ['Phenotips_ID', 'Batch_Received', 'BTRIS Preferred Name', 'Observation Value']]
    df.rename(columns={'BTRIS Preferred Name' : 'Test',
                       'Observation Value': 'CC_Value'}, inplace=True)

    # only include ones that have been sequenced
    df.dropna(inplace = True)
    df = df[df['Batch_Received'].str.contains("BATCH")]

    # clean up some of the test names
    df['Test'] = df['Test'].str.replace("^HLA-", "", regex=True)
    df['Test'] = df['Test'].str.replace("^HLA ", "", regex=True)    
    df['Test'] = df['Test'].str.replace("\*.*$", "", regex=True)
    df['Test'] = df['Test'].str.replace(",.*$", "", regex=True)
    df['Test'] = df['Test'].str.replace(" .*$", "", regex=True)
    df['Test'] = df['Test'].str.replace("Cw", "C", regex=True)
#    df.Test.unique()
    
    df['CC_Value'] = df['CC_Value'].str.replace("\*", "", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace(" ", "", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace(":[0-9]", "", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace("^1$", "01", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace("^3$", "03", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace("^6$", "06", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace("^7$", "07", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace("^8$", "08", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace("^9$", "09", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^,', '', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace(',,', ',', regex=True)
    df.drop(df[df.CC_Value == ""].index, inplace=True)
    df.drop(df[df.CC_Value == "SEEBELOW"].index, inplace=True)
    df = df[~df.CC_Value.str.contains('HLAALLELESARETESTEDBYPCR')]
    
    df['CC_Value'] = df['CC_Value'].str.replace("^9$", "09", regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace(r'(\d\d)\d{1,4}', r'\1', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace(r'(\d\d)/\d\d', r'\1', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace(r'^(\d)$', r'0\1', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^35,35$', '35', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^11,11$', '11', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^03,03$', '03', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^12,12$', '12', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^15,15$', '15', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^04,04$', '04', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^02,02$', '02', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^05,05$', '05', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^06,06$', '06', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^07,07$', '07', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^24,24$', '24', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^44,44$', '44', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^13,13$', '13', regex=True)
    df['CC_Value'] = df['CC_Value'].str.replace('^Feb.*$', '02', regex=True)
    
    return(df)

# Read in the HLA typing from the CSI batch processing
def read_csi(d, pDF):
    # BATCH15/HLA/P0000858/hla/R1_bestguess_G.txt
    # Create the full file names
    #pDF['hlafile'] = d + pDF['Batch_Received'] + "/HLA/" + pDF['Phenotips_ID'] + "/hla/R1_bestguess_G.txt"
    first = True

    #for f in pDF['Dir'].tolist(): 
    for index, row in pDF.iterrows():
        f = d + "/" + str(row['Batch_Received']) + "/HLA/" + str(row['Phenotips_ID']) + "/hla/R1_bestguess_G.txt"

        if os.path.isfile(f):
            # print("Unable to locate HLA file {}".format(f))
            # next
            # print("Found HLA file {}".format(f))
            df = pd.read_csv(f, header=0, encoding='utf-8', sep='\t')
            df['Phenotips_ID'] = row['Phenotips_ID']

            if first:
                allDF = df
                first = False
            else:
                allDF = pd.concat([allDF, df])

    allDF = allDF[ ['Phenotips_ID', 'Locus', 'Chromosome', 'Allele']]
    allDF['Allele'] = allDF['Allele'].str.replace("^.*\*","", regex=True)
    allDF['Allele'] = allDF['Allele'].str.replace(":.*$","", regex=True)
#    print(allDF.head())
    
    x = allDF.groupby(['Phenotips_ID', 'Locus'])['Allele'].apply(lambda x: ','.join(x.sort_values().unique())).reset_index()
    x.rename(columns={'Locus' : 'Test',
                       'Allele': 'CSI_Value'}, inplace=True)
    print(x.head())
    return(x)
    

####################################
# 
# Main 
#
####################################

def main():

    #
    # Usage statement
    #
    parseStr = 'Compare HLA typing from CC and CSI pipeline\n\n\
    Usage:\n\
        hla_CC_vs_wes.py -i infile \n\n\
    Example:\n\
        hla_CC_vs_wes.py -i Lab78706121.csv \n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--cc_file', required=False, nargs='?', type=argparse.FileType('r'), 
                        default='/Volumes/DIR/software/resources/known_HLA/hladata_PTips.csv', help='Input file containing CC HLA typing')
    parser.add_argument('-o', '--outfile', required=False, action='store', type=str, 
                        default="out.csv", help='Output file')
    parser.add_argument('-d', '--csi_dir', required=False, action='store', type=str, 
                        default="/Volumes/DIR/CIDR_DATA_RENAMED", help='Directory containing CSI batch processing')

    args = parser.parse_args()
    cc_file = args.cc_file
    indir = args.csi_dir
    outfile = args.outfile

    ## Import and clean the CC input file
    ccDF = read_CC(cc_file)
    print(ccDF.head())
    
    
    ## Import, concatenate and clean the CSI pipeline HLA typing 
    csiDF = read_csi(indir, ccDF[['Phenotips_ID', 'Batch_Received']])
    print(csiDF.head())

    ## Merge and compare CC and CSI
#    ccDF = df
    bothDF = pd.merge(ccDF, csiDF, how='inner', on=['Phenotips_ID', 'Test'])
    bothDF['Same'] = (bothDF['CC_Value']==bothDF['CSI_Value']).astype(bool)
    bothDF.drop_duplicates(inplace=True)
    bothDF['Duplicates'] = bothDF.duplicated(subset=['Phenotips_ID', 'Test'], keep=False)
    
#    bothDF.loc[bothDF['Duplicates'] == True][['Phenotips_ID', 'Test', 'CC_Value', 'CSI_Value', 'Same']]
    bothDF.loc[bothDF['Same'] == False][['Phenotips_ID', 'Test', 'CC_Value', 'CSI_Value']]
    
    bothDF['Same'].value_counts()
    bothDF['Test'].value_counts()
    len(bothDF['Phenotips_ID'].unique())
    bothDF['Duplicates'].value_counts()
    
    pd.crosstab(bothDF['Test'], bothDF['Same'])
    
if __name__ == '__main__':
    main()

