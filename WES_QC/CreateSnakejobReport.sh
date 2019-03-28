#!/bin/bash

$1

module load snakemake

snakemake --report $1.html --snakefile Multiqc.snakemake
