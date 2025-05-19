#!/bin/bash

SECONDS=0

# STEP1 Download fastq files

# download fastq files (still need to gunzip and move them)
#wget -P ./fastq_files https://figshare.com/ndownloader/articles/24182520/versions/2
#mkdir -p ./transcriptome_files
#mv ./fastq_files/GCF_000146045.2_R64_rna.fna.gz ./transcriptome_files/

# STEP2 Run fastqc and MultiQC

#mkdir -p ./fastqc_output
#fastqc ./fastq_files/*.gz -t 4 -o ./fastqc_output/

multiqc ./fastqc_output/ -d -o ./multiqc_output



