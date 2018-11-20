#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 5 2018
Updated on:

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

ncs_cnv_ids.py
    Reads the NIAID Centralized Sequencing CNV analysis output files
    and replaces the CIDR Array ID with the Phenotips ID 
    and exports it again

v 1.0 - initial code version.

"""
__author__ = 'Susan Huse'
__date__ = 'October 5, 2018'
__version__ = '1.0'

#import sys
import os
#import re
#import datetime
import pandas as pd
#from yaml import load
import argparse
from argparse import RawTextHelpFormatter
from ncbr_huse import send_update, err_out, pause_for_input
 
####################################
# 
# Functions 
#
####################################
# Test file exists or error out
def test_file(f):
    if not os.path.isfile(f):
        err_out("Error: unable to locate input file:  " + f, log)
    else:
        return(True)

# Create the prompt to be used by get_bsi_data
def bsi_prompt_text(bsifile):
    stars = "******************************************************"
    return('\n'.join(["", stars, stars, "Please log into BSI and create a new report that includes:",
    "CIDR Array ID and PhenotipsId\n",
    "These are in Vial Fields, and Subject Tables subject_ncs, respectively",
    "Use the following IDs as search criteria ('CIDR Array ID' =@(like) *inserted IDs*).",
    'Run the report and save the file as: "{}" in the current directory'.format(bsifile),
    '\nPlease enter "y" to print out the list of IDs to use, or "q" to quit.\n']))

# Get the Phenotips information for the subjects in this batch and previous batches
def get_bsi_data(bsifile, theIDs):
    test_file(bsifile)

    # Prompt and import the first set of BSI data - subjects released in this batch
    pause_for_input(bsi_prompt_text(bsifile), 'y', 'q')
    print("\n" + "\n".join(theIDs))
    pause_for_input('\nPlease enter "y" to continue when you have created {}, or "q" to quit.\n'.format(bsifile), 'y','q')

    bsi = import_bsi_phenotips(bsifile)
    return(bsi)

# Import the CNV file header information
def import_cnv_header(filename):
    header_count = 6

    with open(filename) as f:
        header=f.readlines()[0:header_count]
    return("".join(header))

# Import the CNV file column data (below the header info)
def import_cnv(filename, refid):
    headcnt = 5
    colcnt = 6

    # Read it in, use dtype string to avoid float adding lots of extra digits
    df = pd.read_csv(filename, header=headcnt, sep='\t', index_col=0, dtype=str)
    df.index = df.index.astype('str')

    # Fix the weirdness with the columns
    theColumns = df.columns.values.tolist()[1:]
    df = df.iloc[:, : colcnt]
    df.columns = theColumns

    # Remove the reference DNA, e.g., NA12878
    df = df[df.index.str.contains('-'+refid+'-') == False]
    return(df)

# Import the BSI phenotips information
def import_bsi_phenotips(f):
    test_file(f)

    # import and convert to dataframe, using field names
    phenotips = pd.read_csv(f, header=2)
    
    # drop the row count column
    phenotips = phenotips.drop(columns=['OBS'])

    return(phenotips)
    
# Update the IDs from the CIDR Array ID to the Phenotips ID
def remap_ids(df, bsi):
    # Get the column names and add Phenotips to the beginning
    theColumns = ['PhenotipsId'] + df.columns.values.tolist()

    # Merge the two dataframes
    ## SUE!!  Error check for Array IDs that don't return a PhenotipdsIDS
    bothDF = pd.merge(df, bsi, how='left', left_index=True, right_on='CIDR Array ID')
    # Drop the column we don't want 
    bothDF.drop(columns=['CIDR Array ID'], inplace=True)
    # Reset the column order
    bothDF = bothDF[theColumns]

    return(bothDF)

# Write out the header info and the updated columnar data to outfile
def write_new(outfile, header, df):
    # First write out the header information
    f = open(outfile, "w")
    f.write(header + '\n')
    f.close()

    # now write out the column data
    df.to_csv(outfile, sep='\t', header=True, index=False, mode='a')

####################################
# 
# Main 
#
####################################

def main():
    #
    # Usage statement
    #
    parseStr = 'Reads a CNV results file containing Array IDs, and \n\
    outputs a new version replacing Array IDs with Phenotips IDs for GRIS.\n\n\
    Usage:\n\
        ncs_cnv_ids.py -i infile -o outfile -b bsifile\n\n\
    Example:\n\
        ncs_cnv_ids.py -i myInCNV -o myOutCNV\n'

    
    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', required=True, nargs='?', type=argparse.FileType('r'), default=None,
                        help='Input file containing CNV results and Array IDs')
    parser.add_argument('-o', '--outfile', required=True, action='store', type=str, default=None,
                        help='Output file for important results')
    parser.add_argument('-b', '--bsifile', required=False, action='store', type=str, default="bsi.csv",
                        help='Output file from BSI providing Array IDs and Phenotips IDs, default=bsi.csv')
    parser.add_argument('-r', '--refid', required=False, action='store', type=str, default="NA12878",
                        help='Reference DNA sample included with each batch, default=NA12878')

    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    bsifile = args.bsifile
    refid = args.refid

    # 
    # Import the data
    #
    header_info = import_cnv_header(infile.name)
    df = import_cnv(infile.name, refid)

    #
    # Get the Phenotips information for the subjects released in this batch
    #
    bsi = get_bsi_data(bsifile, list(set(df.index.tolist())))

    #
    # Remap the IDs
    #
    df = remap_ids(df, bsi)

    #
    # Write out the new file
    #
    write_new(outfile, header_info, df)

if __name__ == '__main__':
    main()



