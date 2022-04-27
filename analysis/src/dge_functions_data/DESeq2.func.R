library("DESeq2")

#' Perform filtering and differential expression analysis on the read count matrix
#' 
#' @param dds A DESeqDataSet object, which could be created by function 'DESeqDataSetFromMatrix'.
#' @param contrast A string vector with three elements: the column name of group information, name of case group and control group.
#' @param prefix output file prefix.
#' @param minimum_total_count numeric. Minimum total count required.
#' @return NULL.
#' @note The function outputs text file for DEGs statistics and a pdf file for visualization of relevant genes
#' @examples
#' countData <- matrix(rnbinom(10000,mu=200,size=20),ncol=4)
#' rownames(countData) <- paste0("gene", 1:dim(countData)[1])
#' condition <- factor(c("A","A","B","B"))
#' dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
#' contrast_vector <- c("condition", "B", "A")
#' analWithDESeq2(dds, contrast_vector, "./testOutput")
analWithDESeq2 <- function(dds, contrast, prefix, minimum_total_count=15){
	#relevel to make control as reference group
	#dds$condition <- relevel(dds$condition, ref=case_control_names[2])
  library("apeglm")
	#filter out miRNAs with total read count <= 10
	dds <- dds[ rowSums(counts(dds)) > minimum_total_count, ]

	#find DE miRNAs
	dds <- DESeq(dds)
	disp <- dispersions(dds)
	names(disp) <- rownames(dds)
	#get the analysis result
	#res <- results(dds, contrast=contrast)
	res <- lfcShrink(dds, coef=2, type="apeglm")
	#order the result by adjusted p value
	resOrdered <- res[order(res$padj),]
	
	#output the result to a text file
	write.table(as.data.frame(resOrdered), sep="\t", file=paste(prefix, "DESeq2", contrast[2], "vs", contrast[3], "txt", sep="."), col.names=NA, quote=F)
	write.table(disp, sep="\t", file=paste(prefix, "disp_DESeq2", contrast[2], "vs", contrast[3], "disp.txt", sep="."), col.names=NA, quote=F)
	#set file to save the graphs
	pdf(paste(prefix, "DESeq2", contrast[2], "vs", contrast[3], "pdf", sep="."), width=20, heigh=10)

	#set layout of the file
	par(mfrow=c(3,1))
	#plotMA(res, main="DESeq2", ylim=c(-3,3))

	#resMLE <- lfcShrink(dds, coef="condition_ALS_vs_control", type="apeglm")

	#plotMA(resMLE, main="DESeq2", ylim=c(-3,3))
	#loop through each gene with padj <0.05
	
	for(gene in rownames(resOrdered)[!is.na(resOrdered$padj) & resOrdered$padj<0.05]){
		
		#plot their count, group by their condition
		#plotCounts(dds, gene=gene, intgroup=contrast[1], las=2)

	}
	
	dev.off()
	resOrdered
}
