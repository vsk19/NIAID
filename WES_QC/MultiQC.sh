#!/bin/bash

#ExampleRun: bash MultiQC.sh "/path/to/BATCH" "BATCH_QC" "snakejobs" {Note: Remove Quotes and do not include "/" after first argument}


cd $1
mkdir $2
cd $1/$2
mkdir $3
cd $1/$2

cp /hpcdata/dir/SCRIPTS/NCBR_github/NIAID/WES_QC/Multiqc.snakemake /hpcdata/dir/SCRIPTS/NCBR_github/NIAID/WES_QC/MultiQC.sh /hpcdata/dir/SCRIPTS/NCBR_github/NIAID/WES_QC/run.json /hpcdata/dir/SCRIPTS/NCBR_github/NIAID/WES_QC/cluster.json $1/$2

module load snakemake
CLUSTER_OPTS="qsub -pe threaded {cluster.threads} -l h_vmem={cluster.mem} -l h_rt=24:00:00 -wd $1/$2/$3"
snakemake -j 30 --cluster-config cluster.json --cluster "$CLUSTER_OPTS" --keep-going --snakefile Multiqc.snakemake > log.txt 2>&1 &
