# NCBR WES QA/QC Pipeline


**Whole Exome Sequencing** requires a preliminary quality control check. `wes_qc.snakefile` is a collection of tools put together using Snakemake.

**USAGE**:
```bash
# Note: Remove Quotes and do not include "/" after first
bash WESQC "/path/to/BATCH" \
           "BATCH_QC" \
           "snakejobs"
```


##### Tools include the following:
- FastQC
- QualiMap
- Samtools Flagstats
- GATK Select Variants
- GATK Variant Eval
- SNPeff
- BCFtools
- MultiQC

**If using this pipeline, please acknowledge [NIAID Collaborative Bioinformatics Resource](https://ncbr.ncifcrf.gov/).**
