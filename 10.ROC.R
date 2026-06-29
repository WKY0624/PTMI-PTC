library("glmnet")
library("survival")

inputFile="1.3train_expTime_DSS.txt"     
inputFile="1.3test_expTime_DSS.txt"     
geneFile="AI_1.3Genes.txt"

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
rt$DSS=rt$DSS/365   
gene<-read.table(geneFile, header=T, sep="\t", check.names=F)    
sameGene=intersect(colnames(rt),gene$x)
data=rt[,c('DSS','Event',sameGene)]

trainScore=read.table("trainrisk.txt",header=T, sep="\t", check.names=F, row.names=1)
trainScore=trainScore$RS
trainScore=log2(trainScore+1) 
outCol=c("DSS","Event",sameGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"High","Low"))
outTabtr=cbind(data[,outCol],RiskScore=as.vector(trainScore),Risk=risk)
write.table(cbind(id=rownames(outTabtr),outTabtr),file="trainRisk.txt",sep="\t",quote=F,row.names=T)


testScore=read.table("testrisk.txt",header=T, sep="\t", check.names=F, row.names=1)
testScore=testScore$RS
testScore=log2(testScore+1)   
outCol=c("DSS","Event",sameGene)
risk=as.vector(ifelse(testScore>median(testScore),"High","Low"))
outTabte=cbind(data[,outCol],RiskScore=as.vector(testScore),Risk=risk)
write.table(cbind(id=rownames(outTabte),outTabte),file="testRisk.txt",sep="\t",quote=F,row.names=T)

totalrisk=rbind(outTabtr,outTabte)
write.table(totalrisk,file="totalRisk.txt",sep="\t",quote=F,row.names=T)

rt=read.table('testRisk.txt', header=T, sep="\t", check.names=F)
library(tidyverse)
library(survivalROC)
library(timeROC)
ROC_rt=timeROC(T=rt$DSS, delta=rt$Event,    
               marker=rt$RiskScore, cause=1,
               weighting='aalen',
               times=c(1,2,3,4,5,6,7,8,9,10), ROC=TRUE)
auc_text=c()
you_roc <- survivalROC(Stime=rt$DSS,
                       status = rt$Event,
                       marker = rt$RiskScore,
                       predict.time = 5,   
                       method = "KM")

cutoff_5years <- you_roc$cut.values[which.max(you_roc$TP-you_roc$FP)]
cutoff_5years
y1 <- you_roc$TP[you_roc$cut.values==cutoff_5years]   
x1 <- you_roc$FP[you_roc$cut.values==cutoff_5years]   

pdf(file='ROC_test.pdf',width=5,height=5)
plot(you_roc$FP,you_roc$TP, xlab="", ylab="", col='white')
abline(h = seq(0, 1, by=0.1), v = seq(0, 1, by=0.1), col = gray(0.95))
par(new=T)

lines(ROC_rt$FP[,1], ROC_rt$TP[,1], col='#b4daa8', lwd=2)
lines(ROC_rt$FP[,3], ROC_rt$TP[,3], col='#80acf9', lwd=2)
lines(ROC_rt$FP[,5], ROC_rt$TP[,5], col='#E99c93', lwd=2)

legend('bottomright',  
       c(paste0('1-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[1]), ')'),
         paste0('3-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[3]), ')'),
         paste0('5-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[5]), ')')
       ),
       col=c('#b4daa8','#80acf9','#E99c93'), lwd=2, bty = 'n')  

arrows(x0=0.45, y0=0.35,x1=x1,y1=y1-0.2,
       length = 0.08, angle = -20, code = 2,col = "#E57259", lwd = 1.5, lty = 1)
text(0.45,0.3, labels = paste("Cutoff value: ",round(cutoff_5years,3)), col='#E57259', cex=0.9, font=2)

dev.off()

