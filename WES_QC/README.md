# NCBR WES QA/QC Pipeline


**Whole Exome Sequencing** requires a preliminary quality control check. `wes_qc.snakefile` is a collection of tools put together using Snakemake.

**USAGE**:
```bash
# Note: Remove Quotes and do not include "/" after first

bash WESQC "/path/to/BATCH" "BATCH_QC" "snakejobs"

BATCH_QC = A Directory where ouput from Multiqc.snakemake are stored
snakejobs = A Directory where logs of snakemake run are stored
```

**DEPENDENCIES**:

Python2.7
Pip
Multiqc1.7

Install Python2.7 and Pip to the local directory and using pip install MultiQC1.7

**Installation**:

Enter the following commands to download and extract Python 2.7 to your hosting account.
 
 mkdir ~/python
        
 cd ~/python
        
 wget http://www.python.org/ftp/python/2.7.2/Python-2.7.2.tgz
        
 tar zxfv Python-2.7.2.tgz
        
 find ~/python -type d | xargs chmod 0755

 cd Python-2.7.2

**Once extracted you can use the following commands to configure and install Python**

 ./configure --prefix=$HOME/python
        
 make
        
 make install

Modify the .bashrc

For your local version of python to load you will need to add it to the .bashrc file.

 vim ~/.bashrc
        
 Press i 

 Enter:
        export PATH=$HOME/python/Python-2.7.2/:$PATH
        
Write the changes (press ESC) and close vim:
        :wq
        
Press Enter
        source ~/.bashrc

**Download Pip from**

wget https://bootstrap.pypa.io/get-pip.py

type: python get-pip.py --user (This will install pip to your local directory (.local/bin))

set PATH variable for pip in .bashrc file

PATH=$PATH:~/.local/bin

**Download Multqc using following command**

pip install multiqc

-------End-------
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
