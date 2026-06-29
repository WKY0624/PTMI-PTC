library(survival)
library(survminer)
library(timeROC)


inputFile="trainRisk.txt"  
rt=read.table(inputFile, header=T, sep="\t", check.names=F)

ROC_rt=timeROC(T=rt$DSS, delta=rt$Event,    
               marker=rt$TYMS, cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)

pdf(file="TYMS.pdf", width=5, height=5.5)
plot(ROC_rt,time=1,col=mycol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=mycol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=mycol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',  
       c(paste0('1-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[1]), ')'),
         paste0('3-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[2]), ')'),
         paste0('5-Year (AUC=',sprintf("%.03f",ROC_rt$AUC[3]), ')')),
       col=mycol, lwd=2, bty = 'n')   
dev.off()


TCGA = read.table("MRGsall_exp_TPM.txt", header = T, check.names = F, sep = "\t", row.names = 1)  
TCGA = TCGA["TYMS",]
group=sapply(strsplit(colnames(TCGA),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
TCGA = TCGA[,group==0]  
colnames(TCGA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(TCGA)) 
TCGA = as.data.frame(t(TCGA))
TCGA$TYMS = log2(TCGA$TYMS+1)
TCGA$group = ifelse(TCGA$TYMS > median(TCGA$TYMS),"High","Low")

cli = read.table("06.time_DSS.txt", header=T, sep="\t", check.names=F, row.names = 1)
cli$DSS = cli$DSS/365

samSample = intersect(row.names(TCGA),row.names(cli))
rt = cbind(GPR = TCGA[,"group"],cli)

fit <- survfit(Surv(DSS, Event) ~ rt$GPR, data = rt, 
               type = "kaplan-meier", error = "greenwood",
               conf.type = "plain", na.action = na.exclude) 

diff = survdiff(Surv(DSS, Event) ~ rt$GPR, data = rt, na.action = na.exclude)  
pValue = 1 - pchisq(diff$chisq, df = 1)
if(pValue < 0.001){pValue = paste0("P= ",format(pValue, scientific = TRUE))}else{pValue = paste0("P= ",sprintf("%.03f",pValue))}

HR = (diff$obs[1]/diff$exp[1]) / (diff$obs[2]/diff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
HR = sprintf("%.02f",HR)
up95 = sprintf("%.02f",up95)
low95 = sprintf("%.02f",low95)
print(c(HR, low95,up95))

pdf('KM-TYMS.pdf',width =3.5, height =5, onefile=F)
ggsurvplot(fit, 
           data = rt,   
           conf.int = F,  
           pval = pValue,
           pval.size= 5,
           pval.coord = c(0.2,0.7),
           ggtheme = theme_test(base_size = 12, base_line_size = 0.2,base_rect_size = 0.8) +
             theme(plot.title=element_text(hjust=0.5)),  
           title = "TYMS for TCFi", 
           legend.title = 'TYMS',  
           legend.labs=c("High","Low"),  
           palette = mycol, 
           legend = c(0.7, 0.3),  
           font.legend = 12,
           font.main = c(12, "plain","black"),
           font.x = c(12,"plain", "black"),
           font.y = c(12,"plain", "black"),
           font.tickslab = c(12, "plain", 'black'),
           xlab="Time(years)",
           ylab="TCFi probability",
           ylim = c(0.5,1),
           break.time.by = 2,  
           risk.table= T,   
           risk.table.y.text.col = T,    
           risk.table.y.text = T, 
           risk.table.height = 0.2,
           fontsize = 3.5,   
           cumevents=T,
           cumevents.height = 0.2, 
           tables.theme = theme_test(base_size = 8, base_line_size = 0.2, base_rect_size = 0.4)+
             theme(axis.text.y = element_text(size = 10, face = 'bold'),
                   axis.text.x = element_text(size = 10, color = 'black'),
                   legend.title = element_text(size = 10),
                   legend.text = element_text(size = 10),
                   axis.title.x = element_text(size = 10),
                   axis.title.y = element_text(size = 10),
                   axis.ticks = element_line(size = 0.3))
)

dev.off()  

