args <- commandArgs(TRUE)
dir <- args[1]
tissue1 <- args[2]
tissue2 <- args[3]
output_prefix <- args[4]

#dir <- "/Users/suzheng/Documents/suzheng/UNSW/UNSWTasks/2021/velocityRNA/results/SRA_samples/DGE_analysis/SRP032798/figures_tables.ALS_C9.Control.out/"
#tissue1 <- "ALS_C9"
#tissue2 <- "Control"
#outdir <- 

read_one_type_for_deseq2 <- function(feature_type, dir, tissue1, tissue2){
    prefix <- paste0(dir, "/", "DESeq2.", feature_type)
    deseq2_file <- paste(prefix, "DESeq2", tissue1, "vs", tissue2, "txt", sep=".")
    deseq2_out <- read.csv(deseq2_file, header=T, sep="\t", row.names=1)
    colnames(deseq2_out) <- c("logCPM_ori", "logFC", "lfcSE", "stat", "PValue", "FDR")
    deseq2_out$logCPM <- log2(deseq2_out$logCPM_ori+1)
    deseq2_out
}
prefix <- paste0(dir, "/", "EdgeR.", "norm_intron_counts")
edger_file <- paste(prefix, "DEGs", tissue1, tissue2, ".txt", sep="_")
intron_res_edger <- read.csv(edger_file, header=T, sep="\t", row.names=1)

prefix <- paste0(dir, "/", "EdgeR.", "exon_counts_subsampled")
edger_file <- paste(prefix, "DEGs", tissue1, tissue2, ".txt", sep="_")
exon_res_edger <- read.csv(edger_file, header=T, sep="\t", row.names=1)

intron_res_deseq2 <- read_one_type_for_deseq2(feature_type="norm_intron_counts", dir, tissue1, tissue2)
exon_res_deseq2 <- read_one_type_for_deseq2(feature_type="exon_counts_subsampled", dir, tissue1, tissue2)
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
get_stat <- function(intron_res, exon_res, out_prefix="fc.logCPM5"){
    logCPM_thres <- 5
    shared_nonNA_DEG_tables <- get_shared_nonNA_DEG_table(intron_res, exon_res, logCPM_thres)
    intron_edger_shared <- shared_nonNA_DEG_tables$intron_nonNA_shared
    exon_edger_shared <- shared_nonNA_DEG_tables$exon_nonNA_shared
    dge_gene_indices <- (abs(exon_edger_shared$logFC) >0.5 |  abs(intron_edger_shared$logFC) > 0.5) & 
    (exon_edger_shared$PValue < 0.05 |  intron_edger_shared$PValue < 0.05)
    #fc_min <- min(c(exon_edger_shared$logFC[dge_gene_indices], intron_edger_shared$logFC[dge_gene_indices]))
    #fc_max <- max(c(exon_edger_shared$logFC[dge_gene_indices], intron_edger_shared$logFC[dge_gene_indices]))
    e_fc <- exon_edger_shared$logFC[dge_gene_indices]
    i_fc <- intron_edger_shared$logFC[dge_gene_indices]
    #plot(intron_edger_shared$logFC[dge_gene_indices] ~ exon_edger_shared$logFC[dge_gene_indices],
    #  xlab="Log2 fold change in EXONs (FC filtered)", 
    #  ylab="Log2 fold change in INTRONs (FC filtered)",
    #  xlim=c(fc_min,fc_max),
    #  ylim=c(fc_min,fc_max)
    # )
    #lines(c(-10,10), c(-10,10), lty=2, lwd=2)
    #lm_model <- lm(intron_edger_shared$logFC[dge_gene_indices] ~ exon_edger_shared$logFC[dge_gene_indices] )
    #abline(lm_model, col="blue")
    png(paste0(out_prefix, ".abs_fc.png"))
    plot(density(abs(e_fc)-abs(i_fc)))
    dev.off()
    mean_val <- mean(abs(e_fc)-abs(i_fc))
    wt <- wilcox.test(abs(e_fc),abs(i_fc))
    
    list(mean=mean_val, 
         wt_p=wt$p.value, 
         dat=data.frame(e_fc=e_fc, i_fc=i_fc)
        )
}

outprefix <- paste0(output_prefix,".fc_stat.", tissue1, ".", tissue2)
edger_stat <- get_stat(intron_res_edger, exon_res_edger, paste0(outprefix, ".edger"))
write.table(t(c(edger_stat$mean, edger_stat$wt_p)), file=paste0(outprefix, ".edger.stat.txt"), col.names=NA)
write.table(edger_stat$dat, file=paste0(outprefix, ".edger.data.txt"))

deseq2_stat <- get_stat(intron_res_deseq2, exon_res_deseq2, paste0(outprefix, ".deseq2"))

write.table(t(c(deseq2_stat$mean, deseq2_stat$wt_p)), file=paste0(outprefix, ".deseq2.stat.txt"), col.names=NA)
write.table(deseq2_stat$dat, file=paste0(outprefix, ".deseq2.data.txt"))
