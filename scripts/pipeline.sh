#!/bin/bash

SECONDS=0

# STEP0 Set up directory stucture

mkdir -p results/{differential_expression,salmon_mapping,qc/{fastqc_out,multiqc_out/{fastq,salmon_mapping}},visualization}
mkdir -p data/{reference,fastq_raw}
mkdir -p scripts


# STEP1 Download fastq files

wget -P ./data/fastq_raw https://figshare.com/ndownloader/articles/24182520/versions/2
unzip ./data/fastq_raw/2 -d ./data/fastq_raw/ && rm ./data/fastq_raw/2
mv ./data/fastq_raw/GCF_000146045.2_R64_rna.fna.gz ./data/reference && gunzip ./data/reference/GCF_000146045.2_R64_rna.fna.gz


# STEP2 Run fastqc and MultiQC

fastqc ./data/fastq_raw/*.gz -t 4 -o ./results/qc/fastqc_out/
multiqc ./results/qc_control/fastqc_out -d -o ./results/qc_control/multiqc_out/fastq

# STEP3 Run Salmon to map reads

# First, make the transcriptome index
salmon index -t ./data/reference/GCF_000146045.2_R64_rna.fna -i ./data/reference/salmon_transcriptome_index --threads 4

# # Now map the reads
samples=("aerobic_r1" "aerobic_r2" "anaerobic_r1" "anaerobic_r2")

for sample in "${samples[@]}"; do
	salmon quant -i ./data/reference/salmon_transcriptome_index \
		-l A \
		-1 data/fastq_raw/${sample}_1.fq.gz \
		-2 data/fastq_raw/${sample}_2.fq.gz \
		-o results/salmon_mapping/${sample} \
		--threads 4 \
		--validateMappings 
done

# STEP4 Run MultiQC on the mapping files

multiqc -d results/salmon_mapping \
	--outdir results/qc/multiqc_out/salmon_mapping \
	--title "Salmon Mapping QC" \
	--filename salmon_mapping_multiqc_report





