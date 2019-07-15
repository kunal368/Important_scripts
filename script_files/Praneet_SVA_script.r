library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)

args = commandArgs(trailingOnly=TRUE)

path <- args[1]

file1 <- args[2]

file2 <- args[3]

setwd(path)
setwd('/Users/pra7mx/Desktop/Scotts_rna-seq/Human_RNA-SEQ_mapped_to_hg19/Combine_Scott_with_KyleLoh')
file1 <- 'Scotts_and_KyleLoh_RNA-SEQ_TPM_Matrix.txt'
file2 <- 'Sample_Info_Combat.txt'

data_exp <- as.matrix(read.table(file1,sep="\t",header=T,row.names=1, check.names=F))
data_pheno <- as.matrix(read.table(file2,sep="\t",header=T,row.names=1, check.names=F))

pheno <- data.frame(data_pheno)
edata <- data_exp

###### COMBAT ####
mod = model.matrix(~as.factor(type), data=pheno)
mod0 = model.matrix(~1,data=pheno)
n.sv = num.sv(edata,mod,method="leek")
n.sv
svobj = sva(edata,mod,mod0,n.sv=n.sv)
batch = pheno$batch
batch
modcombat = model.matrix(~1, data=pheno)
modcombat
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=FALSE, prior.plots=TRUE)
file <- strsplit(file1,".txt")
outfile <- paste0(file, "_Combat.txt")
write.table(combat_edata,outfile,sep="\t",quote=F)
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")


#pValues = f.pvalue(edata,mod,mod0)
#qValues = p.adjust(pValues,method="BH")

