#!/bin/bash

SECONDS=0

# STEP0 Set up directory stucture

mkdir -p data/{reference,fastq_raw}
mkdir -p mapped_reads_salmon analysis_plots scripts qc/{fastqc_out,multiqc_out/{fastqc,mapped_reads_salmon}} 


## STEP1 Download fastq files

 wget -P ./data/fastq_raw https://figshare.com/ndownloader/articles/24182520/versions/2
 unzip ./data/fastq_raw/2 -d ./data/fastq_raw/ && rm ./data/fastq_raw/2
 rm ./data/fastq_raw/GCF_000146045.2_R64_rna.fna.gz 


## STEP2 Run fastqc and MultiQC

 fastqc ./data/fastq_raw/*.gz -t 4 -o qc/fastqc_out/
 multiqc ./qc/fastqc_out -d -o ./qc/multiqc_out/fastqc

## STEP3 Run Salmon to map reads

## First, get the cDNA references and make the transcriptome index
wget -P data/reference/ https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-61/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz

salmon index -t ./data/reference/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz -i ./data/reference/salmon_transcriptome_index --threads 4

##  Now map the reads

samples=("aerobic_r1" "aerobic_r2" "anaerobic_r1" "anaerobic_r2")

for sample in "${samples[@]}"; do
	salmon quant -i ./data/reference/salmon_transcriptome_index \
		-l A \
		-1 data/fastq_raw/${sample}_1.fq.gz \
		-2 data/fastq_raw/${sample}_2.fq.gz \
		-o mapped_reads_salmon/${sample} \
		--threads 4 \
		--validateMappings 
done

## STEP4 Run MultiQC on the mapping files

 multiqc -d mapped_reads_salmon \
 	--outdir qc/multiqc_out/mapped_reads_salmon \
 	--title "mapped_reads_salmon QC" \
 	--filename mapped_reads_salmon_multiqc_report


# STEP5 Run Analysis_Notebook.qmd to finish analysis. See Notebook for details.
wget -P ./data/reference/ http://sgd-archive.yeastgenome.org/curation/literature/go_slim_mapping.tab

quarto render ../Analysis_Notebook.qmd








