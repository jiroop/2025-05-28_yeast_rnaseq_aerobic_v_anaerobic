## Yeast RNA-seq Analysis: Aerobic vs Anaerobic Growth Conditions

This project analyzes RNA-seq data comparing yeast growth under aerobic and anaerobic conditions using a comprehensive bioinformatics pipeline. The entire analysis pipeline can be run by execution of the bash script ./scripts/pipeline.sh. An analysis notebook describing all steps in the pipeline, the R code used, and a discussion of the analysis done is found in Analysis_Notebook. An html viewable page is [here](https://jiroop.github.io/2025-05-28_yeast_rnaseq_aerobic_v_anaerobic/Analysis_Notebook.html).

Project Overview
The pipeline performs differential gene expression analysis between yeast grown in aerobic versus anaerobic conditions, including:

Quality control of raw sequencing reads
Read mapping and quantification
Differential expression analysis
Data visualization and interpretation

### Pipeline Workflow

1. Data Download: Downloads sequencing reads and yeast transcriptome reference
2. Quality Control: FastQC analysis of raw reads
3. QC Summary: MultiQC report generation for FastQC results
4. Read Mapping: Salmon quasi-mapping for transcript quantification
5. Mapping QC: MultiQC report for Salmon mapping statistics
6. Statistical Analysis: DESeq2 normalization and differential expression analysis in R



### Installation

1. Clone this repository

2. cd 250518_yeast_anaerobic_v_aerobic_expression

3. Create and activate the conda environment:

4. conda env create -f environment.yml

5. conda activate 250518_RNA-seq_yeast_anaerobic_v_aerobic

### Usage

Run the complete pipeline:
./scripts/pipeline.sh

All steps in the pipeline and analysis are described in ./Analysis_Notebook.html (html viewable page [here](https://jiroop.github.io/2025-05-28_yeast_rnaseq_aerobic_v_anaerobic/Analysis_Notebook.html))
