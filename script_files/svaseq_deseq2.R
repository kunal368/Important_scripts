###Upper Quartile Normalization
library(EDASeq)
read.table(file = "Genes_Counts_Matrix_merged.txt", sep="\t")#, quote = FALSE, row.names = FALSE)
counts  <- as.matrix(countData_pFG_DMOG_D6_vs_cond4)

##Filteration and exploratory analysis
filter <- apply(counts, 1, function(x) length(x[x>5])>=2)
counts <- counts[filter,]

#set <- newSeqExpressionSet(as.matrix(counts), phenoData = sampleTable)
UQ_mat <- betweenLaneNormalization(counts, which="upper")
#logTransformed.UQ_mat <- log2(UQ_mat+1)

## Set null and alternative models (ignore batch)
mod1 = model.matrix(~condition)
mod0 = cbind(mod1[,1])
svseq_uq <- svaseq(UQ_mat,mod1,mod0)$sv


########################################Plotting before anf after Upper Quartile Normalization##############
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
#outfile <- paste0("Relative", "_Expression", "_NoNormalization", ".pdf")
png("Relative_Expression_NoNormalization.png",height=1000, width=900)
plotRLE(counts, outline=FALSE, ylim=c(-4, 4), col=colors[x])
dev.off()

png("PCA_before_UQ_normalization.png",height=1000, width=900)
plotPCA(counts, col=colors[x], cex=1.2)
dev.off()

png("Relative_Expression_UQNormalization.png",height=1000, width=900)
plotRLE(UQ_mat, outline=FALSE, ylim=c(-4, 4), col=colors[x])
dev.off()

png("PCA_after_UQ_normalization.png",height=1000, width=900)
plotPCA(UQ_mat, col=colors[x], cex=1.2)
dev.off()
#plotPCA(UQ_mat, col=colors[x], cex=1.2)

##Running deseq on upper-quantile matrix
sampleTable <- data.frame(condition = as.factor(c(rep("Ctl",3),rep("Trt",3))))
rownames(sampleTable) <- colnames(counts)
sampleTable$svseq_uq <- svseq_uq

dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=sampleTable,
                              design=~svseq_uq+condition)

dds <- DESeq(dds)

#pca
vsd <- vst(dds, blind=FALSE)
png("PCA_b4_svaseq.png",height=800, width=900 )
plotPCA(vsd)
dev.off()

result <- as.data.frame(results(dds, contrast=c("condition", "Trt", "Ctl"), alpha = 0.05))
#result <- subset(result, !is.na(result$padj) & result$padj < 0.05)
#result <- subset(result, result$log2FoldChange > 1 | result$log2FoldChange < -1)
result <- result[order(result$log2FoldChange, decreasing = TRUE), ]

#file  <- strsplit(file,split = ".txt")[[1]][1]
write.csv(result, file = "deseq_result_after_svaseq_uq.csv" ) #Writing to csv


####without any batch effect incorporation
dds <- DESeqDataSetFromMatrix(countData=counts,
                             colData=sampleTable,
                             design=~condition)

dds <- DESeq(dds)

result <- as.data.frame(results(dds, contrast=c("condition", "Trt", "Ctl"), alpha = 0.05))
result <- result[order(result$log2FoldChange, decreasing = TRUE), ]
write.csv(result, file = "deseq_result_no_batch_adjustemnt.csv" ) #Writing to csv
