library(limma)                
expFile = "M1_diff_1.3MRGsall_Exp.txt"  
cliFile = "06.time_DSS.txt"

rt=read.table(expFile, header=T, sep="\t", check.names=F)  
rt=as.matrix(rt)
rownames(rt)=rt[,1] 
exp=rt[,2:ncol(rt)] 
dimnames=list(rownames(exp), colnames(exp))  
data=matrix(data=as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)   
group=sapply(strsplit(group,""), "[", 1)    
group=gsub("2", "1", group)   
data=data[,group==0]  
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data)) 
data=avereps(data)   
data=data[rowMeans(data)>0,]  
data=log2(data+1)   
data=t(data)

cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)  

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="1.3total_expTime_DSS.txt",sep="\t",row.names=F,quote=F)   

total = out
trainFile <- read.table("trainSet.txt", header = T, sep = "\t", row.names = 1, check.names = F)
testFile <- read.table("testSet.txt", header = T, sep = "\t", row.names = 1, check.names = F)

train <- total[row.names(trainFile),]
write.table(train,"1.3train_expTime_DSS.txt", sep = "\t", quote = F, row.names = F, col.names = T)

test <- total[row.names(testFile),]
write.table(test,"1.3test_expTime_DSS.txt", sep = "\t", quote = F, row.names = F, col.names = T)


