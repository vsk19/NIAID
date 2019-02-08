#!/bin/bash

#ExampleRun: bash WESQC "/path/to/BATCH" "BATCH_QC" "snakejobs" {Note: Remove Quotes and do not include "/" after first argument}


cd $1
mkdir $2
cd $1/$2
mkdir $3
cd $1/$2

cp /hpcdata/dir/SCRIPTS/MultiQC/Multiqc.snakemake /hpcdata/dir/SCRIPTS/MultiQC/MultiQC.sh /hpcdata/dir/SCRIPTS/MultiQC/run.json /hpcdata/dir/SCRIPTS/MultiQC/cluster.json $1/$2

module load snakemake
CLUSTER_OPTS="qsub -pe threaded {cluster.threads} -l h_vmem={cluster.mem} -l h_rt=24:00:00 -wd $1/$2/$3"
snakemake -j 30 --cluster-config cluster.json --cluster "$CLUSTER_OPTS" --keep-going --snakefile Multiqc.snakemake > log.txt 2>&1 &
