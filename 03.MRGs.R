DESeq=read.table("03.DESeq2.diff.txt", header=T, sep="\t", check.names=F)
DESeq=as.matrix(DESeq)
rownames(DESeq)=DESeq[,1]  
exp=DESeq[,2:ncol(DESeq)]  
dimnames=list(rownames(exp), colnames(exp)) 
DESeq=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames) 

gene=read.table("MRGsall_exp_TPM.txt", header=F, sep="\t", check.names=F)  
sameGene=intersect(as.vector(gene[,1]), rownames(DESeq))  
geneDiff=DESeq[sameGene,]    

out=rbind(ID=colnames(geneDiff),geneDiff)    
write.table(out,file="M1_all_MRGsall_stat.xls",sep="\t",quote=F,col.names=F)

geneDiff2 = as.data.frame(geneDiff)
logFCfilter = 0.379 
fdrFilter = 0.05    
geneFilter=geneDiff2[(abs(as.numeric(as.vector(geneDiff2$log2FoldChange)))>logFCfilter & as.numeric(as.vector(geneDiff2$padj))<fdrFilter),]  

outFilter=rbind(ID=colnames(geneFilter),geneFilter) 
write.table(outFilter, file="M1_diff_1.3MRGsall_stat.xls", sep="\t", row.names=T, col.names=F,quote=F)

Exp <- read.table("05.TPM100.txt", header=T, sep="\t", check.names=F)
Exp=as.matrix(Exp)
rownames(Exp)=Exp[,1] 
exp_Exp=Exp[,2:ncol(Exp)]  
dimnames_Exp=list(rownames(exp_Exp), colnames(exp_Exp))  
Exp=matrix(as.numeric(as.matrix(exp_Exp)), nrow=nrow(exp_Exp), dimnames=dimnames_Exp) 

ExpDiff=Exp[as.vector(rownames(geneFilter)),]
expOut=rbind(ID=colnames(ExpDiff), ExpDiff)
write.table(expOut, file="M1_diff_1.3MRGsall_Exp.txt", sep="\t", col.names=F, quote=F)




