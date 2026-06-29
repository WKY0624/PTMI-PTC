library(limma)
library(sva)
tcgaExpFile="05.TPM100.txt"
rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]    
exp=rt[,2:ncol(rt)]   
dimnames=list(rownames(exp),colnames(exp))   
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)   
tcga=tcga[rowMeans(tcga)>0,]  

tcga6 <- tcga[apply(tcga, MARGIN=1, FUN =  function(xxx) {
  (sum(xxx==0)/length(xxx))<=0.6
           }),]

dim(tcga6)
tcga=log2(tcga+1)

group=sapply(strsplit(colnames(tcga),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tcga=tcga[,group==0]
tcga=t(tcga)
rownames(tcga)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tcga))
tcga=t(avereps(tcga))
tcga[tcga<0]=0

tcgaTab=rbind('Gene ID'=colnames(tcga), tcga)
write.table(tcgaTab, file="total.normalize.txt", sep="\t", quote=F, col.names=F)

total <- tcga
trainFile <- read.table("trainSet.txt", header = T, sep = "\t", row.names = 1, check.names = F)
testFile <- read.table("testSet.txt", header = T, sep = "\t", row.names = 1, check.names = F)

total_train.row <- total[,row.names(trainFile)]
train.nor <- cbind(id=row.names(total_train.row),total_train.row)
write.table(train.nor,"train.normalize.txt", sep = "\t", quote = F, row.names = F, col.names = T)

total_test.row <- total[,row.names(testFile)]
test.nor <- cbind(id=row.names(total_test.row),total_test.row)
write.table(test.nor,"test.normalize.txt", sep = "\t", quote = F, row.names = F, col.names = T)


