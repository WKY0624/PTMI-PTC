library(DESeq2)
library(tidyverse)
library(ggsignif) 
library(RColorBrewer)
library(limma)
library(ggplot2)
library(ggpubr)
library(beepr)
library(gplots)
library(pheatmap)
library(latex2exp)   

rt <- read.table("01.newCounts_477+59.txt",header=T,sep="\t",comment.char="",check.names=F)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp),colnames(exp))
data = matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

data = avereps(data)    
data = data[rowMeans(data)>100,]   
data = data[apply(data, 1, sum) > 0 , ]   
data = as.data.frame(data)
data <- floor(data)  
write.table(data,"01.Counts.txt",sep="\t",quote=F,col.names = NA)


sample.type = sapply(strsplit(colnames(data),"\\-"), "[", 4)
sample.type = sapply(strsplit(sample.type,""), "[", 1)
sample.type = gsub("2", "1", sample.type)
no.normal = length(sample.type[sample.type==1])     
no.tumor = length(sample.type[sample.type==0])     

SampleID = colnames(data)
Type = c(rep(1,no.normal), rep(2,no.tumor))
Group = cbind(SampleID, Type)
Group = as.data.frame(Group)
colnames(Group) = c("id", "Type")
Group$Type = ifelse(Group$Type==1, "Normal", "Tumor")
Group = as.matrix(Group)
rownames(Group) = Group[,1]
Group = Group[,2:ncol(Group)]
Group = as.data.frame(Group)
colnames(Group) = c("condition")
Group$condition <- as.factor(Group$condition)
write.table(Group,"01.Group.txt",sep = "\t", quote = F)


dds <- DESeqDataSetFromMatrix(countData = data,colData = Group,design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Normal")
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
dds <- DESeq(dds)
res <- results(dds)
write.table(res,"03.DESeq2.diff.txt",sep="\t",quote=F,col.names = NA)




vsd <- vst(dds, blind = FALSE)
countData_raw<- assay(dds)
countData_normalized<- assay(vsd)
write.table(countData_normalized,"04.Counts_vst.txt",sep="\t",quote=F,col.names = NA)

n.sample = ncol(data)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))

pdf(file="rawBox.pdf")
boxplot(countData_raw, col = cols, main="expression value",xaxt = "n")
dev.off()
pdf(file="normalBox.pdf")
boxplot(countData_normalized, col = cols,main="expression value",xaxt = "n")
dev.off()
pdf(file="histold.pdf")
hist(countData_raw)
dev.off()
pdf(file="histnew.pdf")
hist(countData_normalized)
dev.off()

count <- read.table('01.Counts.txt', header=T,sep="\t",comment.char="",check.names=F)
tpm <- read.table("01.newTPM_477+59.txt", header=T,sep="\t",comment.char="",check.names=F)
count=as.matrix(count)
rownames(count)=count[,1]
exp=count[,2:ncol(count)]
dimnames=list(rownames(exp),colnames(exp))
count=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
count=as.data.frame(count)

tpm=as.matrix(tpm)
rownames(tpm)=tpm[,1]
exp2=tpm[,2:ncol(tpm)]
dimnames2=list(rownames(exp2),colnames(exp2))
tpm=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)
tpm=as.data.frame(tpm)


sameGene <- intersect(rownames(count), rownames(tpm))
filterTPM <- tpm[sameGene,]
out = rbind(ID = colnames(filterTPM),filterTPM)  
write.table(out,file="05.TPM100.txt",sep="\t",quote=F,col.names=F)



