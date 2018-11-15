#!/usr/bin/env python3

"""
file: RenamePeakBedfile.py
by: Tovah Markowitz
Date: 11/14/18

Purpose: to take a bedfile with duplicated peak names and create a new bed file with
   unique/numbered peak names
   
Inputs: a bed file, any number of columns with at least 4 columns
Outputs: a bed file, same sets of columns as input
"""

##################################################################################
# Modules

import optparse

##################################################################################
# Functions

def readBed(inputName):
    """purpose: to read in bed files and convert into a universal format"""
    # read in input file
    f = open(inputName,'r')
    input = f.readlines()
    f.close()
    # split into a table
    input = [ row.strip().split('\t') for row in input if not row.startswith("#") ]
    input = checkBed(input)
    return(input)

def checkBed(input):
    """purpose: to confirm that the file has the peak identifier in the correct column.
    Also checks columns used for sorting.
    """
    if (len(input[0]) < 4) & (len(input[-1]) < 4):
        print( "File has too few columns to be processed." )
        sys.exit()
    elif not checkNumeric(input[0][1]):
        print( "Not a bed file. Column 2 is not numeric." )
        sys.exit()
    elif checkNumeric(input[0][3]):
        print( "Incorrect file format. Column 4 is numeric." )
        sys.exit()
    elif (input[0][3] == "+") | (input[0][3] == "-") | (input[0][3] == "."):
        print (Incorrect file format. Column 4 contains strand information.")
        sys.exit()
    else:
        return(input)

def checkNumeric(value):
    try:
        int(value)
        return True
    except:
        return False

def Set_Chr_Nr_ (Chr):
    """ Sort by chromosome, found online """
    if Chr: 
        New = Chr[3:]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        elif New == 'M': New = 25
        else: New = int(New)
    else:
        New = 0
    return New

def bedSort(input):
    """sort peaks by chromosome and start position"""
    input = sorted(input,key=lambda x:(Set_Chr_Nr_(x[0]),int(x[1])))
    return(input)

def changePeakName(input):
    """ changes old peak names to incorporate a unique identifier at the end"""
    oldPeakNames = {}
    for peak in input:
        if not peak[3] in list(oldPeakNames.keys()):
            oldPeakNames[peak[3]] = 1
        else:
            oldPeakNames[peak[3]] += 1
        peak[3] = peak[3] + "_" + str(oldPeakNames[peak[3]])
    return(input)

def writeBed(input, outputName):
    """write output files"""
    output = [ '\t'.join(peak) for peak in input ]
    f = open( outputName, 'w' )
    f.write( '\n'.join( i for i in output ) )
    f.close()

################################################################################
# Main

def main():

    desc="""
    This script is designed to take a bed or bed-like file as long as the first four 
    columns are chromosome, start, end, and identifier where the identifier values are 
    not unique. It will take this file, sort by chromosome number and start position, and
    then attach a unique identifier to each current identifier. In other words, the first 
    instance of identifier 'merged_common_peaks' will be renamed 'merged_common_peaks_1',
    the second 'merged_common_peaks_2', etc.
    """
   
    parser = optparse.OptionParser(description=desc)
    # define commandline options
    parser.add_option('-i', dest = 'inputFile', default = '', help = 'The name of the \
input bed file with repetitive peak identifiers.')
    parser.add_option('-o', dest = 'outputFile', default = '', help = 'The name of the \
output bed file with all unique identifiers.')

    (options,args) = parser.parse_args()

    inputName = options.inputFile
    outputName = options.outputFile

    #inputName = "mJP_HCF1_mm_i80_broadPeak_M_above_1.0_biased_peaks.bed"
    input = readBed(inputName)
    input = bedSort(input)
    input = changePeakName(input)
    writeBed(input, outputName)

if __name__ == '__main__':
    main()