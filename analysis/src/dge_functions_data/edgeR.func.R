library(limma)
library(edgeR)
library(RColorBrewer)
#library(mixOmics)
#library(HTSFilter)



#count is the read count matrix
#group is a vector of group of the sample
#prefix the output file prefix
#name of case and control group, will be used as part of output file name
#' Perform filtering, normalization, dispersion estimation and differential expression analysis on the read count matrix
#' 
#' @param counts The numeric matrix of read counts.
#' @param group A factor vector of sample group, only support two groups.
#' @param prefix output file prefix.
#' @param case_control_names A string vector with two elements: name of case group and control group.
#' @param min_count_for_at_least_some_samples numeric. Minimum count required for at least some samples.
#' @param minimum_total_count numeric. Minimum total count required.
#' @return NULL.
#' @note The function outputs 'DEGs' text file for DEGs statistics, 'DEGs' pdf file for data visualization and Design text for input design matrix
#' @examples
#' y <- matrix(rnbinom(10000,mu=200,size=20),ncol=4)
#' group=as.factor(rep(1:2,each=2))
#' anal_counts(y, group, c("A", "B"))

anal_counts <- function(counts, group, prefix, case_control_names, min_count_for_at_least_some_samples=10, minimum_total_count=15){

#file name of DEGs output file, pdf file, design information file
output_DEG <- paste(prefix, "DEGs", case_control_names[1], case_control_names[2], ".txt", sep="_")
output_pdf <- paste(prefix, "DEGs", case_control_names[1], case_control_names[2], ".pdf", sep="_")
output_design <- paste(prefix, "Design", case_control_names[1], case_control_names[2], ".txt", sep="_")
output_disp <- paste(prefix, "disp", case_control_names[1], case_control_names[2], ".disp.txt", sep="_")

#read the data into a DGE object
dge <- DGEList(counts=counts, group=group, remove.zeros = TRUE)

#create the design matrix
design <- model.matrix(~group)

#set the column names of the design matrix
colnames(design) <- sub("group", "", colnames(design))

#output the design matrix
write.table(design, file=output_design)

#index of RNAs to keep, based on the filter
keep <- filterByExpr(dge, design, min.count = min_count_for_at_least_some_samples, min.total.count = minimum_total_count)

#filter the RNAs
dge2 <- dge[keep,,keep.lib.sizes=FALSE]

#calculate normalization factor, to normalized library size
dge2 <- calcNormFactors(dge2)

#estimate the relationship between mean read count and standard variant
dge2 <- estimateDisp(dge2, design=design)

disp <- dge2$tagwise.dispersion
names(disp) <- rownames(dge2)

fit <- glmQLFit(dge2,design)
dgeTest <- glmQLFTest(fit,coef=2)


#do the test to find DEGs
#dgeTest <- exactTest(dge4, pair=case_control_names[c(2,1)])

#get results of all RNAs, ordered by p value
res <- topTags(dgeTest, p.value=1, n = nrow(dgeTest$table))

#output the results
write.table(res, file=output_DEG, sep="\t", quote=F, col.names=NA)
write.table(disp, file=output_disp, sep="\t", quote=F, col.names=NA)

pdf(output_pdf)
oriPar <- par()
par(mfrow = c(4, 3))
pseudoCounts <- log2(dge$counts+1)
for(i in colnames(counts)){
	hist(pseudoCounts[,i], main=i, xlab="log2(count + 1)")
}

par(mfrow = c(1, 1))
par(mar = c(8,4,2,2))
boxplot(log2(dge$counts+1), col="gray", las = 3, cex.names = 1, main="Non-filtered", ylab="log2(count + 1)")
boxplot(log2(dge2$counts+1), col="gray", las = 3, cex.names = 1, main="Filtered", ylab="log2(count + 1)")


cpm2 <- cpm(dge2)
pseudoCounts2 <- cpm(dge2, log = TRUE, prior.count = 1)
boxplot(pseudoCounts2, col = "gray", las = 3, cex.names = 1, main="Non-normalized", ylab="log2(CPM)")




case_i <- group==case_control_names[1]
control_i <- group==case_control_names[2]



mvalues <- rowMeans(pseudoCounts[, case_i]) - rowMeans(pseudoCounts[, control_i])
avalues <- rowMeans(pseudoCounts)

plot(avalues, mvalues, xlab="A", ylab="M", pch=19)
abline(h=0, col="red")

plotMDS(pseudoCounts)

sampleDists <- as.matrix(dist(t(pseudoCounts)))
#par(mar = c(8,4,4,8))
#cim(sampleDists, title="Sample distance matrix")





hist(dgeTest$table[,"PValue"], breaks=50, xlab="P-value")
hist(res$table[,"FDR"], breaks=50, xlab="FDR")

############ definition of positive findings#######
positive <-  res$table$FDR <= 0.05 & abs(res$table$logFC) >= 1

plotSmear(dgeTest, de.tags = rownames(res$table)[which(positive)], cex=0.7)

volcanoData <- cbind(res$table$logFC, -log10(res$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
DEGs4volc <- positive
point.col <- ifelse(DEGs4volc, "red", "black")

plot(volcanoData, pch = 16, col = point.col, cex = 0.5)


y <- cpm(dge, log=TRUE, prior.count = 1)
selY <- y[rownames(res$table)[positive],]
#par(mar = c(8,4,4,8))
#cim(t(selY), title="DEG expression matrix")

dev.off()
res
}


