
####### Reading Command Line Arguments #######
args = commandArgs(trailingOnly=TRUE)

TPM_file_name <- args[1]
NumberofBatchA <- as.numeric(args[2])
NumberofBatchB <- as.numeric(args[3])

TPM_merged_mat <- as.matrix(read.table(TPM_file_name,sep="\t",header=TRUE,row.names=1,check.names=F))

###pca before batch effect removal
data <- TPM_merged_mat[apply(TPM_merged_mat, 1, function(x) !all(x==0)),] #Removing Genes/Entities with zero variance
logTransformed.dat = log2(data+1)
#data.t <- t(data)
#samples <- factor(colnames(logTransformed.dat))
pca <- prcomp(t(logTransformed.dat), center=T, scale. = T)


library(factoextra)
fviz_pca_ind(pca, axes = c(1,2),
             col.ind = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,# Avoid text overlapping
             title = "PCA log-tranformed TPM: Before batch removal")

####Batch removal
batch <- c(rep("A",NumberofBatchA),rep("B",NumberofBatchB))
dat_batch_removed  <- limma::removeBatchEffect(logTransformed.dat, batch)

#pca after batch removal
pca2 <- prcomp(t(dat_batch_removed), center=T, scale. = T)


library(factoextra)
fviz_pca_ind(pca2, axes = c(1,2),
             col.ind = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,# Avoid text overlapping
             title = "PCA log-tranformed TPM: After batch removal")