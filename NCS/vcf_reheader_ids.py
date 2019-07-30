#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 15:54:54 2019
Updated on

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

vcf_reheader_ids.py
    reads a vcf.gz and a masterkey file and renames each of the samples

v 1.0 - initial code version.
"""
__author__ = 'Susan Huse'
__date__ = 'April 10, 2019'
__version__ = '1.0'
__copyright__ = 'No copyright protection, can be used freely'

import gzip
import os
import sys
import re
import datetime
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from ncbr_huse import send_update, err_out, pause_for_input

####################################
# 
# Define Functions
#
####################################
# Test file exists or error out
def test_file(f):
    if not os.path.isfile(f):
        err_out("Error: unable to locate input file:  " + f, log)
    else:
        return(True)

# Import CIDR to Phenotips IDs from  the masterkey_batchX.txt file
def import_masterkey(f):
    test_file(f)

    # import and convert to series
    samplekey = pd.read_csv(f, header=0, sep='\t', encoding='utf-8')
    samplekey = samplekey.set_index('CIDR_Exome_ID')
    #samplekey = samplekey['Phenotips_ID']

    return(samplekey)

######################################
#   
# Functions for writing output files
#
######################################


# Output the three additional files
def write_files(filenames, familyDF, receivedIDs, seqDF, unseqDF):
    #filenames 0: masterkey, 1: batchinfo, 2: seqrped
    receivedDF = familyDF.loc[familyDF['CRIS_Order#'].isin(receivedIDs)]

    #seqFamDF = familyDF[familyDF['Phenotips_ID'].isin(seqIDs)]
    seqDF = pd.concat([receivedDF, seqDF])

    # Write masterkey
    master_fields = ['CRIS_Order#',
                   'Phenotips_ID', 
                   'CIDR_Exome_ID',
                   'Batch_Sent']
    seqDF[master_fields].to_csv(filenames[0], sep='\t', header=True, index=False)
    
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
    parseStr = 'Reads the CIDR Exome IDs in the input VCF file and creates a sampleID file to reheader with\n\
    bcftools reheader --samples <your file> -o <output> <your input file>\n\n\
    Usage:\n\
        vcf_reheader_ids.py -b batch \n\n\
    Example:\n\
        vcf_reheader_ids.py -b 7\n\
        vcf_reheader_ids.py -b 7 -p\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-k', '--masterkey', required=False, nargs='?', type=argparse.FileType('r'), default=None, 
                help='masterkey file, e.g., masterkey_batch12.txt')
    parser.add_argument('-v', '--vcf', required=False, nargs='?', type=argparse.FileType('r'), default=None, 
                help='gzipped VCF file to be reheadered')

    args = parser.parse_args()
    masterkey = args.masterkey
    vcf = args.vcf

    #####################################
    #
    # Set up the variables and the log file
    #
    #####################################

    ## Set up the log file
    #thedate = str(datetime.datetime.now()).split()[0]
    #thedate = re.sub("-","",thedate)
    #global log 
    #log = open('csi_to_gris' + '.log', 'a')
    #log.write('\n' + str(datetime.datetime.now()) + '\n')
    #log.write(' '.join(sys.argv) + '\n')
    #log.write('vcf_reheader_ids.py version ' + __version__ + '\n\n')
    #log.flush()
        
    #####################################
    ##
    ## Read the masterkey file and the vcf line
    ##
    #####################################

    samplekey = import_masterkey(masterkey_file)

    ## SUE!!
    # read in the VCF, 
    # grab the line "^#CHROM", split #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  125771-0224060968
    # split, intersect with the values in the samplekey dataframe, 
    # df.reindex(["Z", "C", "A"]) on the list of values
    # print out the phenotips IDs

    ######################################
    #   
    # Close out and clean up
    #
    ######################################
    send_update("\nvcf_reheader_ids.py successfully completed", log)
    send_update(str(datetime.datetime.now()) + '\n', log)
    log.close()

if __name__ == '__main__':
    main()

