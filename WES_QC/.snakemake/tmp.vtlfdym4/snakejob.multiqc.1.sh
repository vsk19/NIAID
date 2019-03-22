#!/bin/sh
# properties = {"type": "single", "rule": "multiqc", "local": false, "input": [], "output": ["NIAID_Report.html"], "wildcards": {"dir": "NIAID"}, "params": {"patterns": "/hpcdata/dir/CIDR_DATA_RENAMED/references/multiqc_config_file.yaml"}, "log": [], "threads": 1, "resources": {}, "jobid": 1, "cluster": {"threads": "4", "mem": "24G", "vmem": "8G"}}
cd /hpcdata/scratch/boddapatia2/NIAID/WES_QC && \
/sysapps/cluster/software/Anaconda3/5.3.0/envs/snakemake-5.4_env/bin/python3.6 \
-m snakemake NIAID_Report.html --snakefile /hpcdata/scratch/boddapatia2/NIAID/WES_QC/Multiqc.snakemake \
--force -j --keep-target-files --keep-remote \
--wait-for-files /hpcdata/scratch/boddapatia2/NIAID/WES_QC/.snakemake/tmp.vtlfdym4 --latency-wait 5 \
 --attempt 10 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --mode 2  --allowed-rules multiqc  && touch "/hpcdata/scratch/boddapatia2/NIAID/WES_QC/.snakemake/tmp.vtlfdym4/1.jobfinished" || (touch "/hpcdata/scratch/boddapatia2/NIAID/WES_QC/.snakemake/tmp.vtlfdym4/1.jobfailed"; exit 1)

