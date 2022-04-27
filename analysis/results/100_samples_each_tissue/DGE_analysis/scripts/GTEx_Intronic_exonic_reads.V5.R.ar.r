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
plotc = function(x1,x2,
                               ylim=c(min(x2),max(x2)),
                               xlim=c(min(x1),max(x1)),
                               xlab="",ylab="",main="") {
     
    df <- data.frame(x1,x2)
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(x)[1,] + 1L
    cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
    df$col <- cols[df$dens]
    plot(x2~x1, data=df[order(df$dens),], 
         ylim=ylim,xlim=xlim,pch=20,col=col,
         cex=2,xlab=xlab,ylab=ylab,
         main=main)
}

dir <- paste0(out_dir, "/merged_counts/")
setwd(dir)
external_files_folder <- paste0(out_dir, "/scripts/")
id_mapping <- read.csv(paste0(external_files_folder, "/ENSG_Genename_mapping.txt"), sep="\t", header=F)

get_gene_names <- function(ensemblIDs){
    id_mapping[match(ensemblIDs, id_mapping[,10]), 12]
}

read_counts_file <- function(counts_file){
    tran0 <- read.csv(counts_file, header=T, sep="\t", check.names=F)
    #row.names(tran0) <- as.data.frame(strsplit(tran0$Geneid, ":"))[2,]
    row.names(tran0) <- tran0$Geneid
    tran1 <- tran0[,-(1:6)]
    tran_len <- tran0[,c(1,6)]
    list(counts=tran1, tran_len=tran_len)
}

get_trimed_mean <- function(x, quantiles=c(0.1,0.9)){
#	print(length(x))
    cutoffs <- quantile(x, probs = quantiles, na.rm = T)
    mean(x[x>=cutoffs[1] & x<=cutoffs[2]])
}

    set.seed(8)

read_data_of_one_tissue <- function(exon_file, intron_file, normalization=T){
    set.seed(8)
exon_rf <- read_counts_file(exon_file)
intron_rf <- read_counts_file(intron_file)
exon_tran_counts <- exon_rf$counts
intron_tran_counts <- intron_rf$counts
#exon_tran_counts <- exon_tran_counts[,-(1:4)]
#intron_tran_counts <- intron_tran_counts[,-(1:4)]
exon_tran_lens <- exon_rf$tran_len
intron_tran_lens <- intron_rf$tran_len

#Do not remove any low expressed genes, as they may have high expression in other tissues.
high_cnt_exon <- apply(exon_tran_counts, 1, function(x){sum(x>10)>=length(x)*0.0})
high_cnt_intron <- apply(intron_tran_counts, 1, function(x){sum(x>5)>=length(x)*0.0})
high_count_transcripts <- union(rownames(exon_tran_counts)[high_cnt_exon], rownames(intron_tran_counts)[high_cnt_intron])

trans_with_long_exons_introns <- intersect(rownames(exon_tran_lens)[exon_tran_lens[,2] > 200], rownames(intron_tran_lens)[intron_tran_lens[,2] > 200])
ei_len_ratios <- exon_tran_lens[trans_with_long_exons_introns,2]/intron_tran_lens[trans_with_long_exons_introns,2]
trans_with_good_ratios <- trans_with_long_exons_introns[ei_len_ratios >0.001 & ei_len_ratios < 10]

shared_transcripts <- intersect(rownames(exon_tran_counts), rownames(intron_tran_counts))
shared_samples <- intersect(colnames(exon_tran_counts), colnames(intron_tran_counts))
shared_transcripts <- intersect(high_count_transcripts, shared_transcripts)
shared_transcripts <- intersect(trans_with_good_ratios, shared_transcripts)

exon_tran_counts0 <- exon_tran_counts[shared_transcripts, shared_samples]
intron_tran_counts0 <- intron_tran_counts[shared_transcripts, shared_samples]
exon_tran_lens0 <- exon_tran_lens[shared_transcripts,2]
intron_tran_lens0 <- intron_tran_lens[shared_transcripts,2]


#subsample exonic reads, to make exonic read counts roughly equal  to intronic read counts
exon_counts_subsampled <- sapply(1:dim(exon_tran_counts0)[2], function(sample_ind){
    binom_p <- get_trimed_mean(intron_tran_counts0[,sample_ind])/get_trimed_mean(exon_tran_counts0[,sample_ind])
    sapply(exon_tran_counts0[,sample_ind], function(one_expr_val){rbinom(1, one_expr_val,  binom_p)})
}
)
colnames(exon_counts_subsampled) <- colnames(exon_tran_counts0)
rownames(exon_counts_subsampled) <- rownames(exon_tran_counts0)   
if(normalization){
    exon_tran_counts0[exon_tran_counts0==0] <- 0.2
    intron_tran_counts0[intron_tran_counts0==0] <- 0.2
    exon_counts_subsampled[exon_counts_subsampled==0] <- 0.2

    y <- DGEList(counts=exon_tran_counts0)
    y <- calcNormFactors(y, method="TMM")

    exon_tran_counts1 <- t(t(exon_tran_counts0)/y$samples$norm.factors)
    intron_tran_counts1 <- t(t(intron_tran_counts0)/y$samples$norm.factors)
    exon_counts_subsampled <- t(t(exon_counts_subsampled)/y$samples$norm.factors)
}else{
    exon_tran_counts0[exon_tran_counts0==0] <- 0.2
    exon_tran_counts1 <- exon_tran_counts0
    intron_tran_counts1 <- intron_tran_counts0
}

high_conf_cells <- (intron_tran_counts1 + exon_tran_counts1) >100

ie_ratios_raw <- intron_tran_counts1/exon_tran_counts1
#ie_ratios_raw[!high_conf_cells] <- NA

ie_ratios <- ie_ratios_raw * exon_tran_lens0/intron_tran_lens0

#low_NA_transcripts_4_ratios <- rowSums(is.na(ie_ratios)) < dim(ie_ratios)[2]*0.1

#conf_ie_ratios <- ie_ratios[low_NA_transcripts_4_ratios,]

list(norm_exon_counts=exon_tran_counts1,
     norm_intron_counts=intron_tran_counts1,
     ie_ratios_norm_by_len=ie_ratios,
     exon_counts_subsampled=exon_counts_subsampled
    )
}

dge_out <- list()
tissue_exp_data <- list()


#for(tissue in c("Kidney-Cortex", "Liver", "Muscle", "Brain-Cerebellum")){
#for(tissue in c("Kidney-Cortex", "Muscle")){
for(tissue in c(tissue1, tissue2)){
    tissue_exp_data[[tissue]] <- read_data_of_one_tissue(paste0( tissue, ".exon.merged.txt"), paste0(tissue, ".intron.merged.txt"), normalization=F)
}

source(paste0(external_files_folder,"/edgeR.func.R"))
source(paste0(external_files_folder, "/DESeq2.func.R"))


get_dge <- function(tissue1, tissue2, feature_feature_type){
    
    dat1 <- tissue_exp_data[[tissue1]][[feature_type]]
    dat2 <- tissue_exp_data[[tissue2]][[feature_type]]

    counts <- cbind(dat1, dat2)
    tissues <- factor(c(
        rep(tissue1, dim(dat1)[2]), 
        rep(tissue2, dim(dat2)[2])
    ))
    case_control_names <- c(tissue1, tissue2)
    edger_out <- anal_counts(counts=counts,
                group=tissues, 
                prefix=paste0("EdgeR.", feature_type), 
                case_control_names=case_control_names)

    counts_int <- apply(counts, 2, as.integer)
    #counts_int <- apply(counts, 2, ceiling)
    rownames(counts_int) = rownames(counts)
    dds <- DESeqDataSetFromMatrix(counts_int, DataFrame(tissues), ~ tissues)
    contrast_vector <- c("tissues", tissue1, tissue2)
    deseq2_out <- analWithDESeq2(dds, contrast_vector, paste0("DESeq2.", feature_type))
    colnames(deseq2_out) <- c("logCPM_ori", "logFC", "lfcSE", "stat", "PValue", "FDR")
    deseq2_out$logCPM <- log2(deseq2_out$logCPM_ori+1)
    list(edger_out=edger_out,
         deseq2_out=deseq2_out)
    
}

convert_counts_to_integers <- function(counts, ceiling=F){
        if(ceiling){
            counts_int <- apply(counts, 2, ceiling)
        }else{
            counts_int <- apply(counts, 2, round)
        }
        rownames(counts_int) = rownames(counts)
        counts_int
}


#tissue1 <- "Kidney-Cortex"
#tissue2 <- "Muscle"

#tissue1 <- "Liver"
#tissue2 <- "Brain-Cerebellum"

subDir <- paste("figures_tables", tissue1, tissue2, "out", sep=".")
dir.create(file.path(out_dir, "/DGE_out_ar/"), showWarnings = FALSE)
print(file.path(paste0(out_dir, "/DGE_out_ar/"), subDir))
dir.create(file.path(paste0(out_dir, "/DGE_out_ar/"), subDir), showWarnings = FALSE)
setwd(file.path(paste0(out_dir, "/DGE_out_ar/"), subDir))


subsample_exon_counts <- function(tissue){
    norm_exon_counts <- tissue_exp_data[[tissue]]$norm_exon_counts
    norm_intron_counts <- tissue_exp_data[[tissue]]$norm_intron_counts
    norm_exon_counts_int <- convert_counts_to_integers(norm_exon_counts)

    exon_counts_subsampled <- sapply(1:dim(norm_exon_counts_int)[2], function(sample_ind){
        binom_p <- get_trimed_mean(norm_intron_counts[,sample_ind])/get_trimed_mean(norm_exon_counts_int[,sample_ind])
        sapply(norm_exon_counts_int[,sample_ind], function(one_expr_val){rbinom(1, one_expr_val,  binom_p)})
    }
    )
    colnames(exon_counts_subsampled) <- colnames(norm_exon_counts_int)
    exon_counts_subsampled
}


for(feature_type in c("norm_intron_counts", "norm_exon_counts", "exon_counts_subsampled")){
#for(feature_type in c("norm_intron_counts", "exon_counts_subsampled")){
    cmp_name <- paste(tissue1, tissue2, feature_type, sep=".")
    dge_out[[cmp_name]] <- get_dge(tissue1, tissue2, feature_feature_type)
}

pdf("Comparison_of_DGE_analysis_tools.pdf")
intron_cmp_name <- paste(tissue1, tissue2, "norm_intron_counts", sep=".")
#exon_cmp_name <- paste(tissue1, tissue2, "exon_counts_subsampled", sep=".")
exon_cmp_name <- paste(tissue1, tissue2, "norm_exon_counts", sep=".")



compare_two_tools_for_one_feature <- function(edger_out, deseq2_out, main){
    exon_edger_out <- edger_out
    exon_deseq2_out <- deseq2_out
    shared_genes <- intersect(rownames(exon_edger_out), rownames(exon_deseq2_out))
    fdr_e <- -log10(exon_edger_out[shared_genes, "FDR"])
    fdr_d <- -log10(exon_deseq2_out[shared_genes, "FDR"])
    plotc(fdr_e, fdr_d, xlim=c(0,100), ylim=c(0,100), xlab="FDR from EdgeR (-log10)", ylab="FDR from DESeq2 (-log10)")

    lfc_e <- exon_edger_out[shared_genes, "logFC"]
    lfc_d <- exon_deseq2_out[shared_genes, "logFC"]
    plotc(lfc_e, lfc_d, 
          xlab="Log2 fold change from EdgeR", 
          ylab="Log2 fold change from DESeq2",
          main=main
         )
    
}

exon_edger_out <- dge_out[[exon_cmp_name]]$edger_out$table
exon_deseq2_out <- dge_out[[exon_cmp_name]]$deseq2_out
compare_two_tools_for_one_feature(exon_edger_out, exon_deseq2_out, "Analysis based on exonic reads")

exon_edger_out <- dge_out[[intron_cmp_name]]$edger_out$table
exon_deseq2_out <- dge_out[[intron_cmp_name]]$deseq2_out
compare_two_tools_for_one_feature(exon_edger_out, exon_deseq2_out, "Analysis based on intronic reads")
dev.off()

add_half_of_non_zero_min <- function(x){
    x + min(x[x!=0])/2
}

get_shared_nonNA_DEG_table <- function(intron_res, exon_res, logCPM_thres){
    intron_edger_out <- intron_res
    exon_edger_out <- exon_res
    intron_edger_out <- intron_edger_out[intron_edger_out$logCPM>logCPM_thres, ]
    exon_edger_out <- exon_edger_out[exon_edger_out$logCPM>logCPM_thres, ]

    shared_genes <- intersect(rownames(intron_edger_out), rownames(exon_edger_out))
    non_NA_genes <- !is.na(intron_edger_out[shared_genes, "FDR"]) & !is.na(exon_edger_out[shared_genes, "FDR"])
    shared_genes <- shared_genes[non_NA_genes]
    intron_edger_shared <- intron_edger_out[shared_genes, ]
    exon_edger_shared <- exon_edger_out[shared_genes, ]
    list(intron_nonNA_shared=intron_edger_shared,
         exon_nonNA_shared=exon_edger_shared
        )
}

visualize_cpm_fdr <- function(intron_res, exon_res, logCPM_thres, prefix){
    pdf(paste0(prefix, ".DEGs.pdf"))
    shared_nonNA_DEG_tables <- get_shared_nonNA_DEG_table(intron_res, exon_res, logCPM_thres)
    intron_edger_shared <- shared_nonNA_DEG_tables$intron_nonNA_shared
    exon_edger_shared <- shared_nonNA_DEG_tables$exon_nonNA_shared
    plotc(exon_edger_shared$logCPM, 
          intron_edger_shared$logCPM, 
          xlim=c(-2,15), 
          ylim=c(-2,15),
          xlab="Counts per million reads in EXONs (log2)",
          ylab="Counts per million reads in INTRONs (log2)"
         )
    lines(c(0,15), c(0,15), lty=2)
    plotc(-log10(add_half_of_non_zero_min(exon_edger_shared$FDR)), 
          -log10(add_half_of_non_zero_min(intron_edger_shared$FDR)), 
          xlab="FDR in EXONs (-log10)",
          ylab="FDR in INTRONs (-log10)"
         )
    fc_min <- min(c(exon_edger_shared$logFC, intron_edger_shared$logFC))
    fc_max <- max(c(exon_edger_shared$logFC, intron_edger_shared$logFC))
    plotc(exon_edger_shared$logFC,
      intron_edger_shared$logFC,
      xlab="Log2 fold change in EXONs", 
      ylab="Log2 fold change in INTRONs",
      xlim=c(fc_min,fc_max),
      ylim=c(fc_min,fc_max)
     )
    lines(c(-10,10), c(-10,10), lty=2, lwd=2)
    
    dge_gene_indices <- (abs(exon_edger_shared$logFC) >1 |  abs(intron_edger_shared$logFC) > 1) 
    fc_min <- min(c(exon_edger_shared$logFC[dge_gene_indices], intron_edger_shared$logFC[dge_gene_indices]))
    fc_max <- max(c(exon_edger_shared$logFC[dge_gene_indices], intron_edger_shared$logFC[dge_gene_indices]))
    plotc(exon_edger_shared$logFC[dge_gene_indices],
      intron_edger_shared$logFC[dge_gene_indices],
      xlab="Log2 fold change in EXONs (FC filtered)", 
      ylab="Log2 fold change in INTRONs (FC filtered)",
      xlim=c(fc_min,fc_max),
      ylim=c(fc_min,fc_max)
     )
    lines(c(-10,10), c(-10,10), lty=2, lwd=2)
    lm_model <- lm(intron_edger_shared$logFC[dge_gene_indices] ~ exon_edger_shared$logFC[dge_gene_indices] )
    abline(lm_model, col="blue")
    write.table(lm_model$coefficients[2], file=paste0(prefix, ".logFC_cmp.txt"), col.names=F, sep="\t")
    #summary(lm_model)
    fdr_delta <- (abs(exon_edger_shared$logFC) - abs(intron_edger_shared$logFC)) * (abs(exon_edger_shared$logFC) + abs(intron_edger_shared$logFC))/2
    #summary(fdr_delta)
    write.table(t(summary(fdr_delta)), append=T, file=paste0(prefix, ".logFC_cmp.txt"), col.names=T, sep="\t")
    fdr_delta_data <- data.frame(fdr_delta=fdr_delta)
    g <- ggplot(fdr_delta_data, aes(x=fdr_delta)) + 
        geom_histogram() + 
        scale_y_log10() +
        scale_x_continuous(name="Delta of exon intron logFC")
    #print(g)
    
    
    selected_indices <- -log10(exon_edger_shared$FDR) > 3 & -log10(intron_edger_shared$FDR) <1
exon_dge_selected <- exon_edger_shared[selected_indices, ]
intron_dge_selected <- intron_edger_shared[selected_indices, ]
    if(sum(selected_indices) > 1){
        plot(exon_dge_selected$logFC,
              intron_dge_selected$logFC, 
              xlim=c(-4,2), 
              ylim=c(-4,2),
              main="Genes only diff expr in EXONs",
              xlab="log2 FC in EXONs",
              ylab="log2 FC in INTRONs"
             )
        rect(-1, -100, 1, 100, density=3, lty=2)
        plot(exon_dge_selected$logCPM,
              intron_dge_selected$logCPM, 
              xlim=c(5,11), 
              ylim=c(5,11),
              main="Genes only diff expr in EXONs",
              xlab="Counts per million in EXONs (log2)",
              ylab="Counts per million in INTRONs (log2)"
             )
    }
one_feature_only_genes <- rownames(exon_dge_selected[abs(exon_dge_selected$logFC) > 1, ])
write.table(get_gene_names(one_feature_only_genes), file=paste0(prefix, ".only_DE_in_exons_big_logFC.txt"), sep="\t", quote=F, row.names=F, col.names=F)
    
    selected_indices <- -log10(intron_edger_shared$FDR) > 3 & -log10(exon_edger_shared$FDR) <1
exon_dge_selected <- exon_edger_shared[selected_indices, ]
intron_dge_selected <- intron_edger_shared[selected_indices, ]
    if(sum(selected_indices) > 1){
        plot(exon_dge_selected$logFC,
              intron_dge_selected$logFC, 
              xlim=c(-4,2), 
              ylim=c(-4,2),
              main="Genes only diff expr in INTRONs",
              xlab="log2 FC in EXONs",
              ylab="log2 FC in INTRONs"
             )
        rect(-100, -1, 100, 1, density=3, lty=2)
        plot(exon_dge_selected$logCPM,
              intron_dge_selected$logCPM, 
              xlim=c(5,11), 
              ylim=c(5,11),
              main="Genes only diff expr in INTRONs",
              xlab="Counts per million in EXONs (log2)",
              ylab="Counts per million in INTRONs (log2)"
             )
}
one_feature_only_genes <- rownames(intron_dge_selected[abs(intron_dge_selected$logFC) > 1, ])
write.table(get_gene_names(one_feature_only_genes), file=paste0(prefix, ".only_DE_in_introns_big_logFC.txt"), sep="\t", quote=F, row.names=F, col.names=F)
dev.off()
    
}
    

visualize_cpm_fdr(dge_out[[intron_cmp_name]]$edger_out$table, dge_out[[exon_cmp_name]]$edger_out$table, 0, "edgeR.logCPM0")
visualize_cpm_fdr(dge_out[[intron_cmp_name]]$edger_out$table, dge_out[[exon_cmp_name]]$edger_out$table, 5, "edgeR.logCPM5")

visualize_cpm_fdr(dge_out[[intron_cmp_name]]$deseq2_out, dge_out[[exon_cmp_name]]$deseq2_out, 0, "DESeq2.logCPM0")
visualize_cpm_fdr(dge_out[[intron_cmp_name]]$deseq2_out, dge_out[[exon_cmp_name]]$deseq2_out, 5, "DESeq2.logCPM5")


tissue1_ieratio_out <- tissue_exp_data[[tissue1]]$ie_ratios_norm_by_len
tissue2_ieratio_out <- tissue_exp_data[[tissue2]]$ie_ratios_norm_by_len
dim(tissue1_ieratio_out)
dim(tissue2_ieratio_out)

shared_nonNA_DEG_tables <- get_shared_nonNA_DEG_table(dge_out[[intron_cmp_name]]$edger_out$table, dge_out[[exon_cmp_name]]$edger_out$table, 5)
tissue1_ieratio_conf <- tissue1_ieratio_out[rownames(shared_nonNA_DEG_tables$exon_nonNA_shared), ]
tissue2_ieratio_conf <- tissue2_ieratio_out[rownames(shared_nonNA_DEG_tables$exon_nonNA_shared), ]
dim(tissue1_ieratio_conf)

tissue1_ieratio_medians <- apply(tissue1_ieratio_conf, 1, median)
tissue2_ieratio_medians <- apply(tissue2_ieratio_conf, 1, median)

pdf("IE_ratio_differences.pdf")
plot(log10(tissue1_ieratio_medians), 
      log10(tissue2_ieratio_medians),
      xlim=c(-6,3),
      ylim=c(-6,3),
      xlab=paste("median i/e ratios in ", tissue1, "(log10)"),
      ylab=paste("median i/e ratios in ", tissue2, "(log10)")
     )
lines(c(-6,3), c(-6,3), lwd=2, lty=2)
outlier_indices <- abs(log10(tissue1_ieratio_medians) - log10(tissue2_ieratio_medians)) > 1
plot(log10(tissue1_ieratio_medians)[outlier_indices], 
      log10(tissue2_ieratio_medians)[outlier_indices],
      xlim=c(-6,3),
      ylim=c(-6,3),
      xlab=paste("median i/e ratios in ", tissue1, "(log10)"),
      ylab=paste("median i/e ratios in ", tissue2, "(log10)")
     )
lines(c(-6,3), c(-6,3), lwd=2, lty=2)
dev.off()

out_data <- cbind(
    log10(tissue1_ieratio_medians), 
    log10(tissue2_ieratio_medians), 
    log10(tissue1_ieratio_medians) - log10(tissue2_ieratio_medians),
    abs(log10(tissue1_ieratio_medians) - log10(tissue2_ieratio_medians)) > 1,
    get_gene_names(names(tissue1_ieratio_medians))
)
colnames(out_data) <- c(tissue1, tissue2, "delta", "delta_gt1", "gene_name")
write.table(out_data, file="IE_ratio_differences.txt", sep="\t", quote=F, row.names=F, col.names=F)


