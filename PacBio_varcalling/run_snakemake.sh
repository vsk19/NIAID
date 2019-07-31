#! /bin/sh

###################
#
# Launching shell script for NIAID PacBio Pipeline
#
###################
module load snakemake/5.5.2
module load python/3.7

##
## Run snakemake
##
echo "Run snakemake"
#mkdir -p snakejobs

if [ "$1" == "npr" ]
then
    snakemake -npr --snakefile Snakefile
fi

if [ "$1" == "process" ]
then
    snakemake --stats snakemake.stats --restart-times 1 --rerun-incomplete -j 100 --cluster-config cluster.json --cluster "sbatch --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname}" --keep-going --snakefile Snakefile > csi_batch_processing.log 2>&1 &
fi
