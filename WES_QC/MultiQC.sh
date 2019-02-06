#!/bin/bash

#ExampleRun: bash WESQC "/path/to/BATCH" "BATCH_QC" "snakejobs" {Note: Remove Quotes and do not include "/" after first argument}


cd $1
mkdir $2
cd $1/$2
mkdir $3
cd $1/$2
cp /hpcdata/scratch/NIAID_scripts/MultiQC/Multiqc.snakemake /hpcdata/scratch/NIAID_scripts/MultiQC/MultiQC.sh /hpcdata/scratch/NIAID_scripts/MultiQC/run.json /hpcdata/scratch/NIAID_scripts/MultiQC/cluster.json $1/$2

module load snakemake
CLUSTER_OPTS="qsub -pe threaded {cluster.threads} -l h_vmem={cluster.mem} -wd $1/$2/$3 walltime=24:00:00"
snakemake -j 30 --cluster-config cluster.json --cluster "$CLUSTER_OPTS" --keep-going --snakefile wes_qc.snakefile > snakemake.output.txt 2>&1 &
