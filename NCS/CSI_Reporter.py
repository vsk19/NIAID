#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 10:34:26 2019
Updated on

Susan Huse
NIAID Center for Biological Research
Frederick National Laboratory for Cancer Research
Leidos Biomedical

CSI_Reporter.py
    Generates report of NIAID CSI consenting and WES metrics from BSI

v 1.0 - initial code version.

"""
__author__ = 'Susan Huse'
__date__ = 'January 3, 2019'
__version__ = '1.0'
__copyright__ = 'No copyright protection, can be used freely'

import sys
import os
import re
import datetime
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter
from ncbr_bsi import read_conf, send_curl, get_bsi_session, get_bsi_name, bsi_query, return_bsi_info
from ncbr_huse import send_update, err_out, pause_for_input
import dominate
from dominate.tags import *
#from requests.auth import HTTPDigestAuth

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

# Return the query result, using the text and value dictionaries
def write_query(k, print_text=True):
    if valDict[k] is None:
        theCount = "Unknown"
    else:
        theCount = str(valDict[k])

    ktext = "{}: {}".format(descrDict[k], theCount)

    if print_text:
        print(ktext)

    return(ktext)

# Return the number of unique Phenotips ID values
def get_count(df, theField = 'Phenotips_ID'):
    return(len(df[theField].unique()))
        
# Write out the html file
def write_html(outfile):
    doc = dominate.document(title='NIAID CSI Report')

    with doc.head:
        #link(rel='stylesheet', href='ncs.css')
        style("""\
            h1 {
                color: slateblue;
                text-align: center;
            }
            
            h2 {
                color: slateblue;
                text-align: center;
            }
            
            table {
                empty-cells: show;
                width: 75%;
            }
            
            th, td {
                padding: 15px;
                text-align: left;
                align-self: center;
            }
            
            tr:nth-child(even) {
                background-color: #f2f2f2;
            }
        """)

    with doc:
       with div(align="center"):
        attr(cls='body')

        today = datetime.datetime.today()
        h1('NIAID Centralized Sequencing Initiative') 
        h2(today.strftime('Patient Throughput Reporting: %B %d, %Y'))

        with table(border='1').add(tbody()):
            therow = tr()
            with therow:
                td(b("Metric"))
                td(b("Patient Count"))
                td(b("Percent of Consented"))

            for k in descrDict:
                if valDict[k] == 0:
                    theCount = "Unknown"
                else:
                    theCount = str(valDict[k])

                therow = tr()
                with therow:
                    td(b(descrDict[k]), br(), defDict[k])
                    td(theCount)
                    td(str(pctDict[k]) + '%')

    f = open(outfile, "w")
    f.write(doc.render())
    f.close()

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
    parseStr = 'Creates an HTML report of NIAID Centralized Sequencing Initiative metrics, \n\
    including numbers of patients consented, sequenced, processed, and finalized.\n\n\
    Usage:\n\
        CSI_Reporter.py \n\n'

    parser = argparse.ArgumentParser(description=parseStr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-s', '--seqr_batch', required=True, type=int, 
                        help='Batch number of last batch uploaded to SEQR')
    parser.add_argument('-d', '--add_date', required=False, type=bool, default=True,
                        help='append current date "_YYMMDD" to the output filename (default=True)')
    parser.add_argument('-o', '--outfile', required=False, action='store', type=str, 
                        default="CSI_Patient_ThruPut.html", 
                        help='Output file for CSI report (default=CSI_Patient_ThruPut_YYYYMMDD.html)')
    args = parser.parse_args()
    outfile = args.outfile
    add_date = args.add_date
    seqr_batch = args.seqr_batch

    #####################################
    #
    # Set up the log file
    #
    #####################################
    # Set up the log file
    thedate = str(datetime.datetime.now()).split()[0]
    thedate = re.sub("-","",thedate)
    global log 

    log = open('CSI_Reporter' + '.log', 'a')
    log.write('\n' + str(datetime.datetime.now()) + '\n')
    log.write(' '.join(sys.argv) + '\n')
    log.write('CSI_Reporter.py version ' + __version__ + '\n\n')
    log.flush()

    ########################################
    #
    # Set the field names, filter values, and dictionaries
    #
    ########################################
    # Query fields of interest
    fields = ['Phenotips ID', 
              'Phenotips Family ID', 
              'CRIS Order #', 
              'CIDR Exome ID', 
              'Batch Sent', 
              'Batch Received', 
              'CRIS Order Status', 
              'Patient Name']

    # Filter by Order status:
    orders_keep = ['Specimen Collected', 'Auto Complete', 'Corrected Results', 'Final Results', 'Hold']
    orders_canceled = ['Canceled by Performing Dept', 'Canceled', 'Cancelled via Patient Discharge']

    # Set up dictionary for writing out results, and the numeric results dictionary
    global descrDict
    descrDict = {
        'total':'Patient in BSI',
        #'named':'Patients in BSI',
        'not canceled':'Patients who have not been canceled',
        'active':'Active patients',
        'sent': 'Patients whose samples have been sent to CIDR for sequencing',
        'received': 'Patients whose samples have been returned from CIDR and logged into BSI',
        'seqr': 'Patients whose samples have been processed and uploaded to SEQR',
        'finalized': 'Patients whose sample results have been finalized'
        }

    global defDict
    defDict = {
        'total':'(Unique Phenotips IDs in BSI)',
        #'named':'(Records with a valid Patient Name)',
        'not canceled':'(CRIS Order Status not Canceled)',
        'active':'(CRIS Order Status = Specimen Collected)',
        'sent': '(Batch Sent)',
        'received': '(Batch Received)',
        'seqr': '(Batch Received now in SEQR)',
        'finalized': '(CRIS Order Status = Final Results)'
        }

    global valDict
    valDict = {}

    global pctDict
    pctDict = {}
        
    ########################################
    #
    # Collection information from BSI, create the output files
    #
    ########################################
    cnf, url_session, url_reports, curl_get = return_bsi_info()
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)
    ## bsi_query(curl_get, url for reports, session ID, fields to display, query value, query field, isequal)

    # Return all records for which the PhenotipsId != "Q" (random filter that will pull all)
    allRecords = bsi_query(curl_get, url_reports, session, fields, "Q", 'Phenotips ID', False)

    ########################################
    #
    # Filter the records and subset for each query
    #
    ########################################
    q = 'total'
    valDict[q] = get_count(allRecords)

    ## Only records that have a name
    #q = 'named'
    #allRecords = allRecords.loc[allRecords['Patient_Name'] != ""]
    #valDict[q] = get_count(allRecords)

    q = 'not canceled'
    allRecords = allRecords[allRecords['CRIS_Order_Status'].isin(orders_keep)]
    valDict[q] = get_count(allRecords)
    #totalText = write_query(q)

    # Remove the Auto Complete from all remaining queries
    allRecords = allRecords[allRecords['CRIS_Order_Status'] != 'Auto Complete']

    q = 'active'
    valDict[q] = get_count(allRecords.loc[allRecords['CRIS_Order_Status'] == 'Specimen Collected'])
    #activeText = write_query(q)

    q = 'sent'
    valDict[q] = get_count(allRecords.loc[allRecords['Batch_Sent'] != ""])
    #sentText = write_query(q)

    q = 'received'
    valDict[q] = get_count(allRecords.loc[allRecords['Batch_Received'] != ""])
    #receivedText = write_query(q)

    q = 'seqr'
    seqrRecords = allRecords.loc[allRecords['Batch_Received'] != ""][['Phenotips_ID', 'Batch_Received']]
    seqrRecords['Batch_Received'] = seqrRecords['Batch_Received'].str.replace("BATCH0?", "", regex=True)
    seqrRecords['Batch_Received'] = pd.to_numeric(seqrRecords['Batch_Received'])
    seqrRecords = seqrRecords.loc[seqrRecords['Batch_Received'] <= seqr_batch ]
    valDict[q] = get_count(seqrRecords)
    #seqrText = write_query(q)

    q = 'finalized'
    valDict[q] = get_count(allRecords.loc[allRecords['CRIS_Order_Status'] == 'Final Results'])
    #finalizedText = write_query(q)

    for k in valDict:
        pctDict[k] = round(valDict[k] / valDict['total'] * 100)

    ######################################
    #   
    # Creating output files (SampleTracking, Masterkey, Pedigree, Batch Info
    #
    ######################################
    # insert the date 
    if add_date:
        fn, fx = os.path.splitext(outfile)
        outfile = fn + "_" + thedate + fx

    write_html(outfile)

    ######################################
    #   
    # Close out and clean up
    #
    ######################################
    send_update("\n{} successfully written by CSI_Reporter.py".format(outfile), log)
    send_update(str(datetime.datetime.now()) + '\n', log)
    log.close()

if __name__ == '__main__':
    main()

