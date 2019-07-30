#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 10:43:58 2018
Updated on

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

csi_cidr_stats.py
    Reads the sample mapping and VCF information and 
    preps the metadata for uploading to GRIS

v 1.0 - initial code version.

Quality Checks:
    *  which individuals have been returned in this batch
    *  each individual belongs to only one family
    *  family IDs based on Phenotips and MRNs map one-to-one
    *  Order IDs are consistent between pedigree, sample mapping, sample key and orders sent
    *  Finds and removes the duplicate and the control entries
    
"""
__author__ = 'Susan Huse'
__date__ = 'September 5, 2018'
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
#from ncbr_huse import run_cmd, run_os_cmd, un_gzip, send_update, err_out
#import csv
#import scipy
#import numpy

from multiqc.utils import config, report
from multiqc.plots import table_object, beeswarm, bargraph

 
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

# Import pedigree information from Excel file
def import_cidr_stats(f, theColumns):
    test_file(f.name)
    
    # import the data
    df = pd.read_csv(f, header=0)
        
    # subset
    df = df.loc[:, theColumns]

    # Remap the first column and set as index
    df.columns.values[0] = "ExomeID"

    df.set_index('ExomeID', inplace=True)
    df = df.to_dict()

    d = {} #  Initialize the new dictionary as an empty dictionary
    for k in df:
        #d[k] = list(df[k].values())
        d[k] = df[k].values()

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
    parseStr = 'Reads an excel spreadsheet of CIDR derived QC data and creates a multiqc-report like swarm plot.\n\n\
    Usage:\n\
        csi_cidr_stats.py -i qc_report_file -o output_html_file \n\n\
    Example:\n\
        csi_cidr_stats.py -i Holland_Release_Set_10_QC_Report.xlsx -o csi_cidr_stats.html\n'

    
    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', required=True, nargs='?', 
                        type=argparse.FileType('r'), default=None, 
                        help='Input CIDR QC Report Excel file, e.g. "Holland_Release_Set_10_QC_Report.xlsx"')
    parser.add_argument('-o', '--outfile', required=False, nargs='?', 
                        type=argparse.FileType('w'), default=None, 
                        help='Output QC Report HTML file, e.g. "Batch10_stats.html"')
    parser.add_argument('-t', '--testmode', required=False, action='store_true', default=False,
                        help='Run in test mode')

    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    testmode=args.testmode

    #####################################
    #
    # Set up the variables and the log file
    #
    #####################################

    # Set up the log file
    thedate = str(datetime.datetime.now()).split()[0]
    thedate = re.sub("-","",thedate)
    global log 
    log = open('csi_cidr_stats' + '.log', 'a')
    log.write('\n' + str(datetime.datetime.now()) + '\n')
    log.write(' '.join(sys.argv) + '\n')
    log.write('csi_cidr_stats.py version ' + __version__ + '\n\n')
    log.flush()

    ####################################
    #
    # Import Excel file
    #
    ####################################
    theColumns = ['SM_TAG',
                  'VERIFYBAM_AVG_DP',
                  'TOTAL_READS',
                  'PCT_PF_READS_ALIGNED_PAIR',
                  'PF_HQ_ERROR_RATE_PAIR', 
                  'PF_HQ_ALIGNED_Q20_BASES_PAIR',
                  'UNMAPPED_READS',
                  'MEAN_TARGET_COVERAGE',
                  'ZERO_CVG_TARGETS_PCT',
                  'PCT_EXC_MAPQ',
                  'PCT_EXC_BASEQ',
                  'PCT_TARGET_BASES_1X',
                  'PCT_TARGET_BASES_2X',
                  'PCT_TARGET_BASES_10X',
                  'PCT_TARGET_BASES_20X',
                  'PCT_TARGET_BASES_30X',
                  'PCT_TARGET_BASES_40X',
                  'PCT_TARGET_BASES_50X',
                  'PCT_TARGET_BASES_100X'
                  ]

    qcDict = import_cidr_stats(infile, theColumns)
    plot = beeswarm.plot(qcDict)
    print(type(plot))

    ######################################
    #   
    # Close out and clean up
    #
    ######################################
    send_update("\ncsi_cidr_stats.py successfully completed", log)
    send_update(str(datetime.datetime.now()) + '\n', log)
    log.close()

if __name__ == '__main__':
    main()



