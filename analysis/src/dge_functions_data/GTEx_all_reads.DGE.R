library(ggfortify)
library(ggplot2)
library(edgeR)
library(repr)
library(RColorBrewer)
library(ggpubr)
library(DT)
library(limma)
library(edgeR)
library(RColorBrewer)
library(DESeq2)
library(session)

args <- commandArgs(T)

#tissue1 <- "ALS_C9"
#tissue2 <- "Control"
#project <- args[1]
tissue1 <- args[1]
tissue2 <- args[2]
out_dir <- args[3]


dir <- paste0(out_dir, "/merged_counts/")
setwd(dir)
external_files_folder <- "analysis/src/dge_functions_data/"
id_mapping <- read.csv(paste0(external_files_folder, "/ENSG_Genename_mapping.txt"), sep="\t", header=F)

get_gene_names <- function(ensemblIDs){
    id_mapping[match(ensemblIDs, id_mapping[,10]), 12]
}

set.seed(8)

dge_out <- list()
tissue_exp_data <- list()

source(paste0(external_files_folder, "/shared_functions.R"))


#for(tissue in c("Kidney-Cortex", "Muscle")){
for(tissue in c(tissue1, tissue2)){
    tissue_exp_data[[tissue]] <- read_data_of_one_tissue(paste0( tissue, ".exon.merged.txt"), paste0(tissue, ".intron.merged.txt"), normalization=F, tissue=tissue)
}


subDir <- paste("figures_tables", tissue1, tissue2, "out", sep=".")
dir.create(file.path(out_dir, "/DGE_out/"), showWarnings = FALSE)
print(file.path(paste0(out_dir, "/DGE_out/"), subDir))
dir.create(file.path(paste0(out_dir, "/DGE_out/"), subDir), showWarnings = FALSE)
setwd(file.path(paste0(out_dir, "/DGE_out/"), subDir))

for(feature_type in c("norm_intron_counts", "norm_exon_counts", "exon_counts_subsampled")){
#for(feature_type in c("norm_intron_counts", "exon_counts_subsampled")){
    cmp_name <- paste(tissue1, tissue2, feature_type, sep=".")
    dge_out[[cmp_name]] <- get_dge(tissue1, tissue2, feature_type)
}

pdf("Comparison_of_DGE_analysis_tools.pdf")
intron_cmp_name <- paste(tissue1, tissue2, "norm_intron_counts", sep=".")
exon_cmp_name <- paste(tissue1, tissue2, "exon_counts_subsampled", sep=".")

exon_edger_out <- dge_out[[exon_cmp_name]]$edger_out$table
exon_deseq2_out <- dge_out[[exon_cmp_name]]$deseq2_out
compare_two_tools_for_one_feature(exon_edger_out, exon_deseq2_out, "Analysis based on exonic reads")

exon_edger_out <- dge_out[[intron_cmp_name]]$edger_out$table
exon_deseq2_out <- dge_out[[intron_cmp_name]]$deseq2_out
compare_two_tools_for_one_feature(exon_edger_out, exon_deseq2_out, "Analysis based on intronic reads")
dev.off()


visualize_cpm_fdr(dge_out[[intron_cmp_name]]$edger_out$table, dge_out[[exon_cmp_name]]$edger_out$table, 0, "edgeR.logCPM0")
visualize_cpm_fdr(dge_out[[intron_cmp_name]]$edger_out$table, dge_out[[exon_cmp_name]]$edger_out$table, 5, "edgeR.logCPM5")

visualize_cpm_fdr(dge_out[[intron_cmp_name]]$deseq2_out, dge_out[[exon_cmp_name]]$deseq2_out, 0, "DESeq2.logCPM0")
visualize_cpm_fdr(dge_out[[intron_cmp_name]]$deseq2_out, dge_out[[exon_cmp_name]]$deseq2_out, 5, "DESeq2.logCPM5")

exon_cmp_name <- paste(tissue1, tissue2, "norm_exon_counts", sep=".")
visualize_cpm_fdr(dge_out[[intron_cmp_name]]$edger_out$table, dge_out[[exon_cmp_name]]$edger_out$table, 5, "AllReads.edgeRCPM5")
visualize_cpm_fdr(dge_out[[intron_cmp_name]]$deseq2_out, dge_out[[exon_cmp_name]]$deseq2_out, 5, "AllReads.DESeq2CPM5")





