library(survival)
library(survminer)
library(timeROC)
library(pROC)
library(glmnet)

rm(list = ls()) 
 

TCGA = read.table("MRGsall_exp_TPM.txt", header = T, check.names = F, sep = "\t", row.names = 1)   
TCGA = TCGA["MGAT4B",]
TCGA = as.data.frame(t(TCGA))
TCGA$MGAT4B = log2(TCGA$MGAT4B+1)

TCGA$group <- ifelse(substr(rownames(TCGA), nchar(rownames(TCGA)) - 2, nchar(rownames(TCGA)) - 2) == "1", "Normal", "Tumor")
GPL96 = read.table("bindGEO_GPL96_remove.txt", header = T, check.names = F, sep = "\t", row.names = 1)
GPL96 = GPL96["MGAT4B",]
GPL96 = as.data.frame(t(GPL96))
GPL96.group = read.csv("group.csv",row.names = 1)
GPL96.merge <- merge(GPL96, GPL96.group, by = "row.names", all = TRUE)
rownames(GPL96.merge) = GPL96.merge[,1]
GPL96.merge = GPL96.merge[,-1]


GPL570 = read.table("bindGEO_GPL570_remove.txt",header = T, check.names = F, sep = "\t", row.names = 1)
GPL570 = GPL570["MGAT4B",]
GPL570 = as.data.frame(t(GPL570))
GPL570.group = read.csv("group570.csv",row.names = 1)
GPL570.merge <- merge(GPL570, GPL570.group, by = "row.names", all = TRUE)
rownames(GPL570.merge) = GPL570.merge[,1]
GPL570.merge = GPL570.merge[,-1]


pdf('MGAT4B.pdf', height = 4, width = 4)

roc1 <- roc(TCGA$group,TCGA$MGAT4B, ci = T) 
roc2 <- roc(GPL96.merge$group,GPL96.merge$MGAT4B, ci = T) 
roc3 <- roc(GPL570.merge$group,GPL570.merge$MGAT4B, ci = T) 

plot(roc1, title=FALSE, lwd=3, legacy.axes = F,
     xlab="1-Specificity", ylab="Sensitivity", col= mycol[1],
     cex.main=1, cex.lab=1, cex.axis=1, font=1,
    grid= F,   
)
plot.roc(roc2, add=TRUE, 
         col = mycol[2], lwd =3)  
plot.roc(roc3, add=TRUE,  
         col = mycol[3], lwd =3)  

legend('bottomright',
       c(paste0('TCGA (AUC=',sprintf("%.03f",roc1$auc[1]),')'),
         paste0('GPL96 (AUC=',sprintf("%.03f",roc2$auc[1]),')'),
         paste0('GPL570 (AUC=',sprintf("%.03f",roc3$auc[1]),')')),
       col= mycol,
       lwd = 0.5, cex = 0.8,bty = 'n', 
       x.intersp = 0.5, 
       y.intersp = 1)
dev.off()


library(reportROC)

reportROC(PCR$group,
          PCR$MGAT4B,
          important = "se",
          plot = TRUE)

ROC.info <- reportROC(TCGA$group,
                      TCGA$MGAT4B,
                      important = "se",
                      plot = F)

write.csv(ROC.info, "ROC.info.TCGA(log2).csv")


GPL96.merge2 = GPL96.merge
GPL96.merge2$group = ifelse(GPL96.merge2$group =='PTC',1,0)

fit1 <- glm(group ~ MGAT4B,
            data = GPL96.merge2,
            family = binomial())  

summary(fit1)

GPL96.merge2$prob <- predict(fit1, 
                             newdata=GPL96.merge2, 
                             type="response")