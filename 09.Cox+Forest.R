library(survival)  
coxPfilter=0.1 
inputFile1="1.3train_expTime_DSS.txt"   
inputFile2="1.3test_expTime_DSS.txt" 
inputFile3="1.3total_expTime_DSS.txt" 

rt=read.table(inputFile1, header=T, sep="\t", check.names=F, row.names=1)
rt$DSS=rt$DSS/365   

outTab=data.frame()
sigGenes=c("DSS","Event")   
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(DSS, Event) ~ rt[,i], data = rt)   
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  
  if(coxP<coxPfilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

write.table(outTab,file="1.3train.uniCox.txt",sep="\t",row.names=F,quote=F)  

uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="1.3train.uniSigExp.txt",sep="\t",row.names=F,quote=F)  

library(forestplot) 
library(readxl) 

options(scipen = 0)   
options(digits = 3)   

rs_forest1 <- read.table('1.3train.uniCox.txt', header = T, check.names = F, sep = '\t')
rs_forest1<-data.frame(rs_forest1)
rs_forest1 <- rs_forest1[order(rs_forest1$pvalue, decreasing =F),]
rs_forest1$`95%CI` <- paste0('[', round(rs_forest1$HR.95L,3) ,'-', round(rs_forest1$HR.95H,3), ']')
rs_forest1$P.value <- ifelse(rs_forest1$pvalue<0.001,format(rs_forest1$pvalue, scientific = T, digits = 3),
                             sprintf("%.03f", rs_forest1$pvalue))

label<-cbind(c("Gene", rs_forest1$id), 
             c("HR",round(rs_forest1$HR,3)),
             c("95% CI", rs_forest1$`95%CI`),
             c("P value", rs_forest1$P.value))

forestplot(labeltext = label, 
           mean = c(NA,rs_forest1$HR),
           lower = c(NA,rs_forest1$HR.95L), 
           upper = c(NA,rs_forest1$HR.95H), 
           is.summary=c(F,F,F,F,F
                        ), 
           
           lineheight = unit(5,'mm'),    
           colgap = unit(5,'mm'),        
           line.margin=unit(6, 'mm'),
           graphwidth = unit(0.3,"npc"), 
           hrzl_lines = list("2" = gpar(lty = 1, lwd = 2, col="black")),
  
           
           zero = NA, 
           lwd.zero = 1, 
           grid = structure(1, gp = gpar(col = "black", lty=2, lwd=1)), 
           
           graph.pos = 2, 
           lty.zero = '',
           lty.ci = 'solid',  
           lwd.ci = 1.3,   
           ci.vertices = T,  
           ci.vertices.height = 0.15,   
           fn.ci_norm="fpDrawDiamondCI",  
           boxsize = 0.3,   
           
           clip=c(0,7),  
           lwd.xaxis = 1.5,     
           lwd.xticks = 0.5,
           xticks = c(0,2,4,6),   
           xlab="Hazard ratio", 
           xlog= F, 
           
           col=fpColors(box = 'firebrick1',  
                        lines = "black",  
                        summary = "yellow",   
                        text = "black",   
                        axes = "black",   
                        hrz_lines = "black",  
           ),
           
           txt_gp=fpTxtGp(label=gpar(fontfamily="",cex=1),   
                          ticks=gpar(fontfamily="",cex=0.8),    
                          xlab=gpar(fontfamily="",cex=1),     
                          title=gpar(fontfamily="",cex=1)  
           ), 
)

