##Deseq2 complete pipeline along with heatmap generation

#Note:
#Control samples should be placed before the Treatment samples
#Samples in counts file and TPM file should be IDENTICAL and in same order.
#There is no need for high pre-filtering thresholds as deseq2 by itself uses Cook's cut-off and independent filtering to remove outlier counts and very low variable genes respectively.

#Input Arguments to command line:
#1. Counts file
#2. TPM file
#3. NumberofControl
#4. NumberofTreatment
#5. Counts: minimum count threshold for pre-filtering
#6. NumberofSamples: minimum number of samples to which count threshold should be applied (for pre-filtering) 


#Output
#1. PCA plot of counts after variance stabilization using vst function of deseq2
#2. Boxplot of Cook's distance for samples. It helps to see if there are too many outlier counts.
#3. Hirerachical clustering Heatmap of sample-sample distance (Euclidean).
#4. DE_result.xlsx: file containg 3 sheets of results. 1. unfiltered 2. filtered just at FDR < 0.05 3. FILTERED both at FDR < 0.05 AND LFC > |1|
#5. 2 STATIC TPM Heatmaps (using pheatmap package) of Differentially Expressed Genes (based on 2nd and 3rd sheets of DE_result.xlsx). Results are placed in Heatmaps folder.
#6. 2 Dynamic Interactive TPM Heatmaps (using heatmaply package) of Differentially Expressed Genes (based on 2nd and 3rd sheets of DE_result.xlsx). Results are placed in Heatmaps folder.

####### Reading Command Line Arguments #######
args = commandArgs(trailingOnly=TRUE)

#path <- args[1]

filenameR <- args[1] ###Counts_file_name
filename_TPM <- args[2] ##tpm_FILE_NAME
NumberofControl <- as.numeric(args[3])
NumberofTreatment <- as.numeric(args[4])
Counts <- as.numeric(args[5])
NumberofSamples <- as.numeric(args[6])
##############################################


#setwd(path)
file <- filenameR
countData <- as.matrix(read.table(file,sep="\t",header=TRUE,row.names=1,check.names=F))
print ("This is how your input data looks")
head(countData)

## if Rounding needed
#countData <- round(countData,0)

print ("Filtering and exploratory data analysis")
## Filtering and exploratory data analysis
filter <- apply(countData, 1, function(x) length(x[x>Counts])>=NumberofSamples)
filtered <- countData[filter,]
#genes <- rownames(filtered)
#aa <- length(genes)
Controlrep <- rep("Ctl",NumberofControl)
Treatmentrep <- rep("Trt",NumberofTreatment)
#x <- as.factor(c(Controlrep,Treatmentrep))
####SETTING THE DIRECTORY######
#setwd(path)



#Expermiental design
sampleTable <- data.frame(condition = as.factor(c(Controlrep,Treatmentrep)))
rownames(sampleTable) <- colnames(filtered)

library(DESeq2)
#Setting up the DESeqDataSet Object and run the DESeq pipeline
dds <- DESeqDataSetFromMatrix(countData=filtered,
                              colData=sampleTable,
                              design=~condition)
dds <- DESeq(dds)
#dds

##Plotting PCA
file  <- strsplit(file,split = ".txt")[[1]][1]

vsd <- vst(dds, blind=FALSE)
png(paste(file, "PCA.png",sep ="_"),height=800, width=900)
plotPCA(vsd)
dev.off()

####################################################################################################################################################################################
							#To see if there are too many outliers or if one sample is consistently higher than others #
###################################################################################################################################################################################


####Sample Outlier Inspection Using Cook's Distance
##Cook’s distance is a measure of how much a single sample is influencing the fitted coefficients for a gene, 
#and a large value of Cook’s distance is intended to indicate an outlier count
png(paste(file, "boxplot_cooks_distance.png",sep ="_"),height=800, width=900)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

####################################################################################################################################################################################
													#Heatmap of the sample-to-sample distances#
###################################################################################################################################################################################
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition#, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png(paste(file, "sample_dist_heatmap.png",sep ="_"),height=800, width=900)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


##function for doing condition vs control
#deseq_Result_cal  <- function(dds, condition, control, out_file){
result <- as.data.frame(results(dds, contrast=c("condition", "Trt", "Ctl"), alpha = 0.05))
result <- result[order(result$log2FoldChange, decreasing = TRUE), ]

result_fdr_05 <- subset(result, !is.na(result$padj) & result$padj < 0.05) #filtering out genes with FDR=NA or FDR < 0.05
result_fdr_05_lfc1 <- subset(result_fdr_05, result_fdr_05$log2FoldChange > 1 | result_fdr_05$log2FoldChange < -1) #filtering further with LFC > |1|


##writing result to csv and xlsx
#write.csv(result, file = paste(file, "DE_result_unfiltered.csv", sep= "_")) #Writing to csv the DE genes at FDR < 0.05

#writing to xlsx
library("openxlsx")
wb = createWorkbook() #Creating workbook

#Adding 3 worksheets
addWorksheet(wb, "All_genes")
addWorksheet(wb, "DE_genes_FDR_0.05")
addWorksheet(wb, "DE_genes_FDR_05_lfc1")
#addDataFrame(result, sheet=sheet1, startColumn=1, row.names=FALSE)
#addDataFrame(dataframe2, sheet=sheet, startColumn=10, row.names=FALSE)

#Writing data

writeData(wb,sheet=1, x=result,rowNames=TRUE,keepNA=TRUE)
writeData(wb, sheet=2, x=result_fdr_05,rowNames=TRUE)#startColumn=1, row.names=FALSE)
writeData(wb, sheet=3, x=result_fdr_05_lfc1,rowNames=TRUE)#startColumn=1, row.names=FALSE)

saveWorkbook(wb, file=paste(file, "DE_result.xlsx", sep= "_"), overwrite=TRUE)

####################################################################################################################################################################################
									#Heatmap Generation from TPM for Differentialy expressed genes#
####################################################################################################################################################################################
#Reading TPM file 
TPM_Data <- as.matrix(read.table(filename_TPM,sep="\t",header=TRUE,row.names=1,check.names=F))

##Checking if Heatmaps directory is already present, if not creating one and changing current dir to Heatmaps dir
mainDir <- getwd()
subDir <- "Heatmaps"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

#Subsetting TPM dataframe for Differentially expressed genes (at 1. FDR < 0.05 only 2. FDR <0.05 AND LFC > |1|)
tpm_fdr_05 <- subset(TPM_Data, rownames(TPM_Data) %in% rownames(result_fdr_05)) #1. FDR < 0.05 only
tpm_fdr_05_lfc1 <- subset(TPM_Data, rownames(TPM_Data) %in% rownames(result_fdr_05_lfc1))#2. FDR <0.05 AND LFC > |1|


#filebase_TPM  <- strsplit(filename_TPM,split = ".txt")[[1]][1]
#Creating static heatmaps using pheatmap

png("TPM_FDR_05_static_heatmap.png",height=800, width=900)
pheatmap(tpm_fdr_05, scale = "row", show_rownames = FALSE, cluster_cols = FALSE)#FDR < 0.05 only
dev.off()

png("TPM_FDR_05_LFC1_static_heatmap.png",height=800, width=900)
pheatmap(tpm_fdr_05_lfc1, scale = "row", show_rownames = FALSE, cluster_cols = FALSE)#DGEs at FDR <0.05 AND LFC > |1|
dev.off()

##Creating dynamic heatmaps using heatmaply
library(heatmaply)

heatmaply(tpm_fdr_05,colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), Colv = FALSE,scale="row",
          showticklabels = c(T,T) ,width = 1000, height = 900, plot_method = "plotly", file = "TPM_FDR_05_dynamic_heatmap.html")

heatmaply(tpm_fdr_05_lfc1,colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), Colv = FALSE,scale="row",
          showticklabels = c(T,T) ,width = 1000, height = 900, plot_method = "plotly", file = "TPM_FDR_05_LFC1_dynamic_heatmap.html")










