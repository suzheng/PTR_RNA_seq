library(edgeR)

args <- commandArgs(T)
#tissue1 <- args[1]
#tissue2 <- args[2]
input_prefix <- args[1]
id_mapping_file <- args[2]
output_prefix <- args[3]
#input_prefix <- "/Users/suzheng/Documents/suzheng/UNSW/UNSWTasks/2021/velocityRNA/results/other_tissues_gene_levels/data/Kidney-Cortex"
meta_file <- "analysis/results/100_samples_each_tissue/DGE_analysis/sample_attributes.txt"
#id_mapping_file <- "/Users/suzheng/Documents/suzheng/UNSW/UNSWTasks/2021/velocityRNA/results/featureCounts_gene_levels/ENSG_Genename_mapping.txt"
#output_prefix <- "/Users/suzheng/Documents/suzheng/UNSW/UNSWTasks/2021/velocityRNA/results/SRA_samples/DGE_analysis/SRP032798/SRP032798.singleTissue"

out_dir <- dirname(output_prefix)
dir.create(file.path(out_dir), showWarnings = FALSE)
setwd(out_dir)

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

id_mapping <- read.csv(id_mapping_file, sep="\t", header=F)

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
high_cnt_exon <- apply(exon_tran_counts, 1, function(x){sum(x>10)>=length(x)*0.1})
high_cnt_intron <- apply(intron_tran_counts, 1, function(x){sum(x>5)>=length(x)*0.1})
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
ie_ratios_raw[!high_conf_cells] <- NA

ie_ratios <- ie_ratios_raw * exon_tran_lens0/intron_tran_lens0

low_NA_transcripts_4_ratios <- rowSums(is.na(ie_ratios)) < dim(ie_ratios)[2]*0.1

conf_ie_ratios <- ie_ratios[low_NA_transcripts_4_ratios,]

list(norm_exon_counts=exon_tran_counts1[low_NA_transcripts_4_ratios,],
     norm_intron_counts=intron_tran_counts1[low_NA_transcripts_4_ratios,],
     ie_ratios_norm_by_len=conf_ie_ratios,
     exon_counts_subsampled=exon_counts_subsampled[low_NA_transcripts_4_ratios,]
    )
}

tissue_exp_data <- read_data_of_one_tissue(paste0(input_prefix, ".exon.merged.txt"), paste0(input_prefix, ".intron.merged.txt"), normalization=T)

conf_ie_ratios <- tissue_exp_data$ie_ratios_norm_by_len

ratio_means <- log10(colMeans(conf_ie_ratios, na.rm=T))
ratio_medians <- log10(apply(conf_ie_ratios, 2, function(x){median(x, na.rm=T)}))
ratio_sds <- log10(apply(conf_ie_ratios, 2, function(x){sd(x, na.rm=T)}))

#ieratios <- tissue_exp_data$ie_ratios_norm_by_len
#ratio_medians <- apply(ieratios, 1, median)
#ratio_means <- apply(ieratios, 1, mean)
#ratio_sds <- apply(ieratios, 1, sd)


exp <- tissue_exp_data$norm_exon_counts

ratio_means_exp_cors <- apply(exp, 1, function(x){cor(x, ratio_means)})
ratio_medians_exp_cors <- apply(exp, 1, function(x){cor(x, ratio_medians)})
ratio_sds_exp_cors<- apply(exp, 1, function(x){cor(x, ratio_sds)})
ratio_mean_d_sds_exp_cors<- apply(exp, 1, function(x){cor(x, ratio_means/ratio_sds)})

#ratio_median_d_exp_cors <- apply(exp, 1, function(x){cor(x, ratio_medians/ratio_sds)})
#ratio_median_d_exp_cor_pvals <- apply(exp, 1, function(x){cor.test(x, ratio_medians/ratio_sds)$p.value})
#ratio_median_d_exp_cor_pvals_adj <- p.adjust(ratio_median_d_exp_cor_pvals, method="bonferroni")

ratio_medians_exp_cor_pvals <- apply(exp, 1, function(x){cor.test(x, ratio_medians)$p.value})
ratio_medians_exp_cor_pvals_adj <- p.adjust(ratio_medians_exp_cor_pvals, method="bonferroni")

length(ratio_medians_exp_cor_pvals_adj[ratio_medians_exp_cor_pvals_adj<0.01])

cor_data <- cbind(get_gene_names(names(ratio_medians_exp_cors)), ratio_medians_exp_cors, ratio_medians_exp_cor_pvals_adj)

write.table(cor_data, file=paste0(output_prefix, ".ieratio_exp_cor.txt"), sep="\t", quote=F, col.names=NA)

meta_df <- read.csv(meta_file, sep="\t", header=T, stringsAsFactors=T)
rownames(meta_df) <- meta_df$SAMPID

shared_samples <- match(names(ratio_medians), meta_df$SAMPID)

numeric_cols <- unlist(lapply(meta_df, is.numeric))  
meta_df1 <- meta_df[shared_samples, numeric_cols]

factor_cols <- unlist(lapply(meta_df, is.factor))  
meta_df2 <- meta_df[shared_samples, factor_cols[]]

aov_pvals <- vector()
aov_names <- vector()
for(i in 1:dim(meta_df2)[2]){    
  tryCatch({ 
   current_p <- summary(aov(ratio_medians~meta_df2[,i]))[[1]][["Pr(>F)"]][1]
   
   if(is.numeric(current_p)){
       aov_names <- c(aov_names, colnames(meta_df2)[i])
       aov_pvals <- c(aov_pvals, current_p)
   }      #print(i)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
}
names(aov_pvals) <- aov_names
aov_ratio_medians_meta_cor_pvals <- aov_pvals

aov_ratio_medians_meta_cor_pvals_adj <- p.adjust(unlist(aov_ratio_medians_meta_cor_pvals), method="bonferroni")
#aov_ratio_medians_meta_cor_pvals_adj
non_na_counts <- apply(meta_df1, 2, function(x)sum(!is.na(x)))
meta_df1_clean <- meta_df1[,non_na_counts > 50]
                       
ratio_medians_meta_cors <- apply(meta_df1_clean, 2, function(x){cor(ratio_medians, x)})
ratio_medians_meta_cor_pvals <- apply(meta_df1_clean, 2, function(x){cor.test(ratio_medians, x)$p.value})

ratio_medians_meta_cor_pvals_adj <- p.adjust(ratio_medians_meta_cor_pvals, method="bonferroni")

write.table(c(ratio_medians_meta_cor_pvals_adj, aov_ratio_medians_meta_cor_pvals_adj), file=paste0(output_prefix, ".meta_cor.txt"), sep="\t", quote=F, col.names=NA)
write.table(ratio_medians, file=paste0(output_prefix, ".ratio_medians.txt"), sep="\t", quote=F, col.names=NA)
