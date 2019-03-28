# NCBR WES QA/QC Pipeline


**Whole Exome Sequencing** requires a preliminary quality control check. `Multiqc.snakefile` is a collection of tools put together using Snakemake.

**USAGE**:
```bash
# Note: Remove Quotes and do not include "/" after first

-->> bash MultiQC.sh "/path/to/BATCH" "BATCH{putbatchnumberhere}_QC" "snakejobs"

BATCH{putbatchnumberhere}_QC = A Directory where ouput from Multiqc.snakemake are stored
snakejobs = A Directory where logs of snakemake run are stored

After successful generation of OUTPUT files, run the following shell script

-->> bash CreateSnakejobReport.sh BATCH{put batch number here} (e.g. bash CreateSnakejobReport.sh BATCH1)

"IF" not successful, following are the reasons:

1. No Samples in VCF file matching the BAM files
2. No available free nodes on the cluster, eventually terminating the snakejobs abruptly without error
```
**Install Dependencies only if not present in your system**
**DEPENDENCIES**:

1. Python2.7
2. Pip
3. Multiqc1.7

"Install Python2.7 and Pip to the local directory and using pip install MultiQC1.7"

**Installation**:

1. Enter the following commands to download and extract Python 2.7 to your hosting account.
 
       mkdir ~/python
        
       cd ~/python
        
       wget http://www.python.org/ftp/python/2.7.2/Python-2.7.2.tgz
        
       tar zxfv Python-2.7.2.tgz
        
       find ~/python -type d | xargs chmod 0755

       cd Python-2.7.2

2. Once extracted you can use the following commands to configure and install Python

       ./configure --prefix=$HOME/python
        
       make
        
       make install

3. For your local version of python to load you will need to add it to the .bashrc file.

   Modify the .bashrc

       vim ~/.bashrc
        
       Press i 

  Enter:
        
       export PATH=$HOME/python/Python-2.7.2/:$PATH
        
4. Write the changes (press ESC) and close vim:
        
        :wq
        
Press Enter
        
        source ~/.bashrc

**Download Pip from**

wget https://bootstrap.pypa.io/get-pip.py
 
        python get-pip.py --user(This will install pip to your local directory (.local/bin))

5. set PATH variable for pip in .bashrc file

       PATH=$PATH:~/.local/bin

6. Download Multqc using following command

       pip install multiqc

7. Goto $HOME/.local/lib/python2.7/site-packages/multiqc/utils/ and edit snpeff argument in search_patterns.yaml from

snpeff:
    contents: 'SnpEff_version'
    max_filesize: 100000    
         
TO        

snpeff:
    contents: 'SnpEff_version'
    max_filesize: 130000

Note: In your multiqc.snakemake, change the path of Multiqc in rule multiqc to your local Multiqc directory that you installed using above commands.

***-------End-------***

##### Tools include the following:
- FastQC
- QualiMap
- Samtools Flagstats
- CollectVariantCallMetrics
- VCFtools
- GATK Select Variants
- GATK Variant Eval
- SNPeff
- BCFtools
- MultiQC

**If using this pipeline, please acknowledge [NIAID Collaborative Bioinformatics Resource](https://ncbr.ncifcrf.gov/).**
