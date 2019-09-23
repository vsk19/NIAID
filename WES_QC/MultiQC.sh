#!/bin/bash

#ExampleRun: bash MultiQC.sh "/path/to/BATCH" "BATCH_QC" "snakejobs" {Note: Remove Quotes and do not include "/" after first argument}


cd $1
mkdir -p $2
cd $1/$2
mkdir -p $3
cd $1/$2

cp /hpcdata/dir/NCBR/NIAID/WES_QC/Variant_metrics_barplot.R /hpcdata/dir/NCBR/NIAID/WES_QC/vcftools_barplot.R /hpcdata/dir/NCBR/NIAID/WES_QC/CreateSnakejobReport.sh /hpcdata/dir/NCBR/NIAID/WES_QC/Multiqc.snakemake /hpcdata/dir/NCBR/NIAID/WES_QC/MultiQC.sh /hpcdata/dir/NCBR/NIAID/WES_QC/run.json /hpcdata/dir/NCBR/NIAID/WES_QC/cluster.json $1/$2

module load snakemake/5.3.0-Python-3.5.5
CLUSTER_OPTS="qsub -sync y -pe threaded {cluster.threads} -l h_vmem={cluster.mem} -l virtual_free={cluster.vmem} -o $1/$2/$3/ -e $1/$2/$3/ -cwd"
snakemake -k --stats snakemake.stats --rerun-incomplete --restart-times 10 -j 100 --cluster-config cluster.json --cluster-sync "$CLUSTER_OPTS" --keep-going --snakefile Multiqc.snakemake > log.snakemake.txt 2>&1 &
