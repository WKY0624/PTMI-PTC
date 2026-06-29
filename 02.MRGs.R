library(limma)  

rt=read.table("05.TPM100.txt", header=T, sep="\t", check.names=F)#Count
rt=as.matrix(rt)
rownames(rt)=rt[,1] 
exp=rt[,2:ncol(rt)] 
dimnames=list(rownames(exp), colnames(exp))  
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames) 

gene=read.table("17inter.txt", header=F, sep="\t", check.names=F)  #
sameGene=intersect(as.vector(gene[,1]), rownames(data))     
geneExp=data[sameGene,]  

out=rbind(ID=colnames(geneExp),geneExp) 
write.table(out,file="MRGsall_exp_TPM.txt",sep="\t",quote=F,col.names=F)
