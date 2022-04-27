# Overview
This is the workflow to perform the bioinformatics analysis in the paper:

Zheng Su etc. Post-transcriptional regulation supports the homeostatic expression of mature RNA.

# System requirements
## Hardware requirements
The workflow requires a PBS Linux cluster and a standard Linux computer with >=16GB RAM to support the in-memory operations.

### Operating system

	Linux 3.10.0-1160.45.1.el7.x86_64 x86_64
	PBSPro version 19.1.3

### Singularity images
Some analyses were performed in Sigularity containers. Please set environmental variable \$SI to be the path of your Sigularity images folder.

	featurecounts.sif
	genometools.sif
	exon_intron_dge_r_packages.sif
	
### Third-party software dependencies
	singularity/3.6.4
	fastqc/0.11.8
	bedtools/2.27.1
	star/2.7.9a
	samtools/1.13
	java/8u45
	conda/4.10.3
	snakemake/6.5.3
	python/3.9.6
	python/2.7.5
	R/4.1.2

### R packages
apeglm v1.14.0,
DESeq2 v1.32.0,
DT v0.19,
edgeR v3.34.1,
ggfortify v0.4.12,
ggplot2 v3.3.5,
ggpubr v0.4.0,
limma v3.48.3,
RColorBrewer v1.1.2,
repr v1.1.3,
session v1.0.3,
dplyr v1.0.7,
ggstatsplot v0.9.1,
ggthemes v4.2.4,
grid v4.1.2,
hexbin v1.28.2,
purrr v0.3.4,
scales v1.1.1

# Installation guide
### Pull singularity images

	singularity pull library://suzheng/ptr_rna/featurecounts:v1.0
	singularity pull library://suzheng/ptr_rna/genometools:v1.0
	singularity pull library://suzheng/ptr_rna/exon_intron_dge_r_packages:v1.0
### Installation of third-party software
Install dependent third-party software following the instruction of the software.

# Bioinformatics Analysis
## Before started
* Most of the bash commands used to perform the analysis are saved in the cmd.sh files, to run them, we usually have to change the working directory to the folder of used cmd.sh file.
* Set an alias using command `alias l="less"`

## Preparation of annotation and reference index files
Working directory:

	srv/scratch/oateslab/share/data/hg38/Gencode/ 
Download gencode.v38.annotation.gff3.gz from [Gencode database](https://www.gencodegenes.org/human/release_38.html) 
 
 
Generate the human gff files for all exons and introns, as well as 5' end only,  3' end only exons and introns, using commands:
 
 	srv/scratch/oateslab/share/data/hg38/Gencode/cmd.sh
which will generate gff files:

 	gencode.v38.GRCh38.exonsOnly.clean.gff3
 	gencode.v38.GRCh38.intronsOnly.exonsSubtracted.gff3.gz
 	gencode.v38.GRCh38.exonsOnly.clean.5prime.gff3
 	gencode.v38.GRCh38.exonsOnly.clean.3prime.gff3
 	gencode.v38.GRCh38.intronsOnly.exonsSubtracted.5prime.gff3
 	gencode.v38.GRCh38.intronsOnly.exonsSubtracted.3prime.gff3
 
STAR2 reference index can be downloaded from:
 
 	gs://gtex-resources/STAR_genomes/STARv275a_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v34_oh100.tar.gz and placed under srv/scratch/oateslab/share/data/hg38/STAR_ref_index/STARv275a_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v34_oh100/
 	
Collapse transcript annotation for RNA-seq QC:

	srv/scratch/oateslab/share/data/non_human_species/collapsed_genes_gtf/cmd.sh
 
## GTEx dataset analysis

### Data download
Working directory for the GTEx dataset analysis:
`analysis/results/100_samples_each_tissue/`

Download the bam and bai files for samples in `sample_names.txt` to directory `out/SAMPLE_NAME/`, in which SAMPLE_NAME is the sample name.

### Quality control

Commands to perform quality control:

	analysis/results/qc/scripts/cmd.sh
	
Collect rnaseqc quality control results, extract QC-failed samples:

	analysis/results/qc/anal_rnaseqc/cmd.sh

Collect fastq QC summaries, extract failed samples:

	analysis/results/qc/anal_fastqc/cmd.sh

### Quantification and differential expression analysis

Perform quantification:
	
	for s in `cat sample_names.txt`;do echo "sh run_one_sample.sh $s >log";done|sh
	
Perform differential expression analysis using all samples or only selecting 5 samples from each tissue:
	
	analysis/results/100_samples_each_tissue/DGE_analysis/cmd.sh

### Using only boundary reads or reads at 5' end or 3' end exons or introns
Working directory:

	analysis/results/stat_with_boundary_reads
	
Commands to perform quantification:

	analysis/results/stat_with_boundary_reads/cmd.sh

Differential expression analysis using boundary reads, working directory:

	analysis/results/stat_with_boundary_reads/DGE_analysis
Commands to perform differential expression analysis:

	analysis/results/stat_with_boundary_reads/DGE_analysis/cmd.sh

Differential expression analysis uisng 5' end or 3' end exonic or intronic reads, working directory:

	analysis/results/stat_with_boundary_reads/eeDGE_analysis
Commands:

	analysis/results/stat_with_boundary_reads/eeDGE_analysis/cmd.sh

## SRA dataset analysis
### Data download and alignment
URLs for download can be found in

	srv/scratch/oateslab/rawData/2021/SRA_RNA_seq/*/*txt
Commands to peform data download, alignment and quantification for the SRA samples:

	analysis/results/SRA_samples/cmd.sh
 
### Quality control
Quality control of GTEx, SRA and non-human dataset was done together using the same set of commands as shown in the GTEx sample quality control step.

### Quantification and differential expression analysis
Perform differential gene expression analysis for the SRA samples:

	analysis/results/SRA_samples/DGE_analysis/cmd.sh

Perform differential gene expression analysis using only 5' end exons or introns for the SRA samples:

	analysis/results/SRA_samples/5eeDGE_analysis/cmd.sh
Perform differential gene expression analysis using only 3' end exons or introns for the SRA samples:

	analysis/results/SRA_samples/3eeDGE_analysis/cmd.sh
 
Meta data (grouping information) for the differential expression analysis, only quality control passed samples included:

	analysis/results/SRA_samples/meta/*.tsv
	
### Using only exon-exon or exon-intron boundary reads

Working directory to perform analysis using only exon-exon or exon-intron boundary reads:
 
 	analysis/results/boundary_reads_only_SRA/
 
Perform quantification:
 	
 	analysis/results/boundary_reads_only_SRA/cmd.sh

Differential expression analysis:

	analysis/results/boundary_reads_only_SRA/DGE_analysis/cmd.sh
 	

## Non-human dataset analysis
### Download annotation files and reference genome files for non-human species
Commands to perform download:

	srv/scratch/oateslab/share/data/non_human_species/anno/cmd.sh

### Generate annotation files for exons and introns, and STAR reference index files
Commands:

	srv/scratch/oateslab/share/data/non_human_species/scripts/cmd.sh
	
Config files for quantification and differential expression analysis pipeline:

	srv/scratch/oateslab/share/data/non_human_species/config/*
	
Perform data download, alignment and quantification:

	analysis/results/non_human_species/cmd.sh
	
### Quality control
Quality control of GTEx, SRA and non-human dataset was done together using the same set of commands as shown in the GTEx sample quality control step.

### Differential expression analysis
Working directory:
	
	analysis/results/non_human_species/DGE_analysis
Meta data files:

	analysis/results/non_human_species/DGE_analysis/meta_*/*.tsv
	
Perform differential expression analysis:

	analysis/results/non_human_species/DGE_analysis/cmd.sh

## Statistics and result visualization
### Collect analysis results for statistics and result visualization
Working directory:
	
	analysis/results/combined_analysis
Commands to collect analysis results:

	analysis/results/combined_analysis/read_counts_intron_exon/cmd.sh
	analysis/results/combined_analysis/tables_cbind/cmd.sh
	analysis/results/combined_analysis/cmd.sh
 
### Statistics and result visualization
 Jupyter notebook to generate Figures 1a and 2a:

	scatter_plot_examples.ipynb


 Jupyter notebook to generate Figure 1b, 2b, 1c, 2c and Figure S2:

	fc_stat_stat.ipynb
  Jupyter notebook to generate Figures 1d, 1e, 1f and 1g:

 	visualize_tissue_diff.ipynb
 Jupyter notebook to generate Figure 3 and Figure S1b:

	delta_fc_cpm_hk_genes.ipynb
 Jupyter notebook to generate Figure S1a:

	read_counts_exon_intron.ipynb
 Jupyter notebook to generate Figure S3:

	reading_strategy_schemic.ipynb
	
# License
This project is covered under the Apache 2.0 License.
        



 
  