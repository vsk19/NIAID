#! /bin/sh

###################
#
# Launching shell script for NIAID CSI batch processing of WES data
#
###################
module load snakemake

##
## Test commandline arguments
##
if [ $# -ne 1 ]; then
    echo " " 
    echo "Requires a single commandline argument: gris, npr, or process"
    echo " " 
    exit
fi

if [ $1 != "gris" ] && [ $1 != "npr" ] && [ $1 != "process" ] ; then
    echo " " 
    echo "Invalid commandline option: $1"
    echo "Valid commandline options include: gris, npr, or process"
    echo " " 
    exit
fi

##
## Get batch and batch number
##
batchdir=`pwd`
batch=`echo $batchdir | sed -e 's/^.*\///' `
echo "BATCH: $batch"
batchnumber=`echo $batch | sed -e 's/BATCH//' -e 's/^0//' `
echo "Processing Batch $batchnumber"

##
## Find the raw directory
##
raw_root="/hpcdata/dir/CIDR_DATA_RAW"
raw_dir="${raw_root}/${batch}/Released_Data_Batch${batchnumber}/Holland_WES_Release_${batchnumber}"

if ! test -d "rawdata"; then
    if test -d $raw_dir; then
        echo "Linking rawdata subdirectory to $raw_dir"
        ln -s $raw_dir rawdata
    else
        echo "Unable to locate raw data directory $raw_dir"
        echo "Exiting"
        exit
    fi
else
    echo "input directory rawdata already exists"
fi

##
## Make the new output directories
##
for i in BAM VCF QC/TARGET QC/UCSC
do
    if ! test -d $i; then
        echo "Creating output directory: $i"
        mkdir -p $i
    else
        echo "output directory $i already exists"
    fi
done


##
## Run csi_to_gris.py
##
if [ "$1" == "gris" ] 
then
    echo "Running csi_to_gris"
    csi_to_gris.py -b $batchnumber
    exit
fi


##
## Run snakemake
##
echo "Run snakemake"

CLUSTER_OPTS="qsub -pe threaded 8 -l h_vmem=32 -l virtual_free=32 -wd $batchdir"

if [ "$1" == "npr" ]
then
    snakemake -npr --snakefile /hpcdata/dir/CIDR_DATA_RENAMED/csi_batch_processing.snakemake
fi

if [ "$1" == "proces" ]
then
    snakemake -k --stats snakemake.stats --rerun-incomplete --restart-times 10 -j 100  --cluster "$CLUSTER_OPTS" --keep-going --snakefile /hpcdata/dir/CIDR_DATA_RENAMED/csi_batch_processing.snakemake > csi_batch_processing.log 2>&1 &
fi
