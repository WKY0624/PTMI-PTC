library(survival)
library(survminer)

bioSurvival=function(inputFile=NULL, outFile=NULL){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F,row.names=1)
  diff=survdiff(Surv(DSS, Event) ~Risk,data = rt)  
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue= paste0("P = ", format(pValue, scientific = TRUE, digits = 3))
  }else{
    pValue=paste0("P = ",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(DSS, Event) ~ Risk, data = rt)  
  
  HR = (diff$obs[1]/diff$exp[1]) / (diff$obs[2]/diff$exp[2])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
  HR = sprintf("%.02f",HR)
  up95 = sprintf("%.02f",up95)
  low95 = sprintf("%.02f",low95)
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,   
                     conf.int.alpha = 0.1,
                     conf.int.style = 'ribbon',
                     palette = mycol,   
                     pval = F,  
                     pval.size = 4,   
                     pval.coord=c(2,0.25),   
                     pval.method = TRUE,  
                     pval.method.size = 4, 
                     pval.method.coord=c(2,0.35),   
                     ggtheme = theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
                       theme(axis.title.y = element_text(size = 10),
                             axis.title.x = element_text(size = 10),
                             axis.text.y = element_text(size = 8),
                             axis.text.x = element_text(size = 9),
                             legend.title = element_text(size = 10),
                             legend.text = element_text(size = 10),
                             axis.ticks = element_line(size = 0.3)
                       ),
                     legend.title="",
                     legend.labs=c("High-risk", "Low-risk"),   
                     legend = c(0.8, 0.2),  
                     font.legend = 10,
                     xlab="Follow-up time (year)",
                     ylab="TCFi probability", 
                     break.time.by = 1,
                     ylim = c(0.5,1),
                     )
  surPlot$plot <- surPlot$plot + 
    annotate("text",x = 1, y = 0.58, hjust = 0, fontface = 2, label = pValue) 
  
  
  
  pdf(file=outFile,onefile = FALSE,width = 3.5, height = 3.5)
  print(surPlot)
  dev.off()
}


bioSurvival(inputFile="trainRisk.txt", outFile="trainSurv.pdf")
bioSurvival(inputFile="testRisk.txt", outFile="testSurv.pdf")
bioSurvival(inputFile="totalRisk.txt", outFile="totalSurv.pdf")
