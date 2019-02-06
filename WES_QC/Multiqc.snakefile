import os
from os.path import join
BAM_ID, = glob_wildcards("../BAM/" + "{ID}.bam")
VCF_ID, = glob_wildcards("../VCF/" + "{ID}.vcf.gz")

configfile:"run.json"

rule all:
    input:fastqc=expand(join("FastQC/{Sample}_fastqc.html"),Sample=BAM_ID),
    	flagstats=expand(join("Flagstats/{Sample}.flagstats"),Sample=BAM_ID),
    	  qualimap=expand(join("{Sample}", "qualimapReport.html"),Sample = BAM_ID),
    	  varianteval=expand(join("VariantEval/{Sample}"),Sample = BAM_ID),
    	  snpeff= expand(join("SNPeff/{Sample}/{Sample}"),Sample = BAM_ID),
          bcftools=expand(join("BCFStats/{Sample}"),Sample = BAM_ID),
          multiqc=join("BatchQC_Report.html"),
 	  cumulativeqc = join("CumulativeQC_Report.html")

rule fastqc:
	input: join("../BAM","{Sample}.bam")
	output: join("FastQC/{Sample}_fastqc.html")
	params: adapters=config['references']['fastqc_adapters']
	threads: 16
	shell: "module load fastqc;fastqc -o FastQC -f fastq --threads {threads} -f bam --contaminants {params.adapters} {input}"

rule qualimap:
 	input: join("../BAM/","{Sample}.bam")
   	output: txt = join("{Sample}","genome_results.txt"), html = join("{Sample}", "qualimapReport.html")
	threads:16
   	params:regions=config['references']['REGIONS'], dir = "{Sample}"
   	shell: "module load qualimap;unset DISPLAY; qualimap bamqc -bam {input} --java-mem-size=48G -c gd hg19 -ip -outdir {params.dir} -gff {params.regions} -outformat HTML -nt {threads} --skip-duplicated -nw 500 -p NON-STRAND-SPECIFIC"

rule samtools_flagstats:
  	input:bam= join("../BAM/","{Sample}.bam")
  	output:join("Flagstats/{Sample}.flagstats")
  	shell: "module load samtools; samtools flagstat {input} > {output}"

rule Gatk_SelectVariants:
        input:selectvariants = expand(join("../VCF/","{VCF}.vcf.gz"),VCF=VCF_ID)
        output:temp(join("{Sample}.vcf.gz"))
        params:genome=config['references']['GENOME'], Sname = "{Sample}"
        shell: "module load GATK/3.7-0-Java-1.8.0_92;java -Xmx48g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R {params.genome} -o {output} -V {input} --sample_name {params.Sname} --excludeNonVariants"

rule bcftools:
		input:"{Sample}.vcf.gz"
  		output:join("BCFStats/{Sample}")
  		shell: "module load bcftools/1.4.1-goolf-1.7.20; bcftools stats {input} > {output}"

rule varianteval:
  		input:vcf = "{Sample}.vcf.gz"
  		output:join("VariantEval/{Sample}")
   		params:genome=config['references']['GENOME'],vcf=config['references']['DBSNP']
   		threads: 2
   		shell:"module load GATK/3.7-0-Java-1.8.0_92;java -Xmx12g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantEval -R {params.genome} -o {output} --dbsnp {params.vcf} --eval {input.vcf} -nt {threads}"

rule snpeff:
  		input:"{Sample}.vcf.gz"
  		output:vcf=join("SNPeff/{Sample}/{Sample}_exome.vcf"),
	           csv = join("SNPeff/{Sample}/{Sample}"),
	           html = join("SNPeff/{Sample}/{Sample}.html")
  		params:genome=config['references']['SNPEFF_GENOME'],effconfig=config['references']['SNPEFF_CONFIG']
  		shell: "module load java/1.8.0_92; java -Xmx12g -jar /hpcdata/scratch/lackjb/snpEff/snpEff.jar -v -canon -c {params.effconfig} -csvstats {output.csv} -stats {output.html} {params.genome} {input} > {output.vcf}"

rule multiqc:
 		input:expand(join("FastQC/{Sample}_fastqc.html"),Sample=BAM_ID),
 			  expand(join("Flagstats/{Sample}.flagstats"),Sample=BAM_ID),
 			  expand(join("{Sample}", "qualimapReport.html"),Sample = BAM_ID),
 			  expand(join("VariantEval/{Sample}"),Sample = BAM_ID), 
 			  expand(join("SNPeff/{Sample}/{Sample}"),Sample = BAM_ID), 
 			  expand(join("BCFStats/{Sample}"),Sample = BAM_ID)
 		output:out1 = "BatchQC_Report.html", out2 = "CumulativeQC_Report.html"
 		shell:"$HOME/.local/bin/multiqc -f -n {output.out1} .;multiqc -f -n {output.out2} ../"
