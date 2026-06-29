library(survival)
library(survminer)


riskFile ="totalRisk.txt"   

risk=read.table(riskFile, header = T, sep='\t', check.names = F)
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


cliFile = "06.clinical477.txt"
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names = 1)

samSample = intersect(row.names(TCGA),row.names(cli))
rt = cbind(TYMS_ = TCGA[,"group"],cli)
rt = cbind(risk,rt)
rt<-na.omit(rt)

rt$TD1 <- rt$TumorDiameter
rt$TD1 = ifelse(rt$TD1 =='0.1~1.0',1,2)

rt$TD4 <- rt$TumorDiameter
rt$TD4 = ifelse(rt$TD4 =='4.1~max',2,1)

rt$HT <- rt$`CombinedDisease`
rt$HT = ifelse(rt$HT =='Lymphocytic thyroiditis',1,0)

rt$`CombinedDisease2` <- rt$`CombinedDisease`
rt$`CombinedDisease2`=ifelse(rt$`CombinedDisease2`=='None',0,
                              ifelse(rt$CombinedDisease2=='Lymphocytic thyroiditis',2,1)) 
rt$TNM2 <- rt$TNM
rt$TNM2 <- ifelse(rt$TNM2 == 'Stage I',12,
                  ifelse(rt$TNM2 =='Stage II',12,34)) 

rt$TNM3 <- rt$TNM
rt$TNM3 <- ifelse(rt$TNM3 == 'Stage I',1,
                  ifelse(rt$TNM3 =='Stage II',2,3)) 

rt$ETE2 <- rt$ETE
rt$ETE2 <- ifelse(rt$ETE2 == 'None',"None","ETE") 

ageless55 <- rt[rt$Age=='<55',c('DSS','Event','TYMS_')] 
agemore55 <- rt[rt$Age=='>=55',c('DSS','Event','TYMS_')]

female <- rt[rt$Gender=='Female',c('DSS','Event','TYMS_')]
male <- rt[rt$Gender=='Male',c('DSS','Event','TYMS_')]

T4 <- rt[rt$`Tcategory`=='T4',c('DSS','Event','Risk')]

N0 <- rt[rt$Ncategory=='N0',c('DSS','Event','TYMS_')]
N1 <- rt[rt$Ncategory=='N1',c('DSS','Event','TYMS_')]
M1 <- rt[rt$`Mcategory`=='M1',c('DSS','Event','Risk')]

TD4_less <- rt[rt$TD4=='1', c('DSS','Event','TYMS_')]
TD4_more <- rt[rt$TD4=='2', c('DSS','Event','TYMS_')]

TD1_less <- rt[rt$TD1=='1', c('DSS','Event','TYMS_')]
TD1_more <- rt[rt$TD1=='2', c('DSS','Event','TYMS_')]

ETE0 <- rt[rt$ETE=='None',c('DSS','Event','TYMS_')]
ETE1 <- rt[rt$ETE2=='ETE',c('DSS','Event','TYMS_')]

histo1 <- rt[rt$`HistologicalType`=='Classical',c('DSS','Event','TYMS_')]
histo2 <- rt[rt$`HistologicalType`=='Follicular variant',c('DSS','Event','TYMS_')]
histo3 <- rt[rt$`HistologicalType`=='Aggressive variants',c('DSS','Event','TYMS_')]

CD0 <- rt[rt$`CombinedDisease2`=='0',c('DSS','Event','Risk')]
CD1 <- rt[rt$`CombinedDisease2`=='1',c('DSS','Event','Risk')]
CD2 <- rt[rt$`CombinedDisease2`=='2',c('DSS','Event','Risk')]

HT0 <- rt[rt$HT=='0',c('DSS','Event','Risk')]

Focal0 <-rt[rt$Multifocality=='Unifocal',c('DSS','Event','TYMS_')]
Focal1 <-rt[rt$Multifocality=='Multifocal',c('DSS','Event','TYMS_')]

BRAF0 <- rt[rt$BRAF=='Wild-type',c('DSS','Event','TYMS_')]
BRAF1 <- rt[rt$BRAF=='Mutation',c('DSS','Event','TYMS_')]

TERT0 <- rt[rt$TERT=='Wild-type',c('DSS','Event','TYMS_')]
TERT1 <- rt[rt$TERT=='Mutation',c('DSS','Event','TYMS_')]

RAS0 <- rt[rt$RAS=='Wild-type',c('DSS','Event','TYMS_')]
RAS1 <- rt[rt$RAS=='Mutation',c('DSS','Event','TYMS_')]

GF0 <- rt[rt$GeneFusion=='Absence',c('DSS','Event','TYMS_')]
GF1 <- rt[rt$GeneFusion=='Presence',c('DSS','Event','TYMS_')]

TNM12 <- rt[rt$TNM2 =='12', c('DSS','Event','TYMS_')]
TNM34 <- rt[rt$TNM2 =='34', c('DSS','Event','TYMS_')]

TNM1 <- rt[rt$TNM3 =='1', c('DSS','Event','TYMS_')]
TNM2 <- rt[rt$TNM3 =='2', c('DSS','Event','TYMS_')]
TNM3_34 <- rt[rt$TNM3 =='3', c('DSS','Event','TYMS_')]

ATA1 <- rt[rt$ATA=='Low', c('DSS','Event','TYMS_')]
ATA2 <- rt[rt$ATA=='Intermediate', c('DSS','Event','TYMS_')]
ATA3 <- rt[rt$ATA=='High', c('DSS','Event','TYMS_')]


diff=survdiff(Surv(DSS, Event) ~ATA,data = rt)  
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue= paste0("P = ", format(pValue, scientific = TRUE, digits = 3))
}else{
  pValue=paste0("P = ",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(DSS, Event) ~ ATA, data = rt)  

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
                   legend.labs=c('High','Intermediate','Low'),  
                   legend = c(0.8, 0.2),   
                   font.legend = 10,
                
                   xlab="Follow-up time (year)",
                   ylab="TCFi",  
                   break.time.by = 1,
                   ylim = c(0.5,1),
)
surPlot$plot <- surPlot$plot + 
  annotate("text",x = 0.2, y = 0.60, hjust = 0, fontface = 1, label = paste0("HR = ",HR,"[",low95,"-",up95,"]")) +
  annotate("text",x = 0.2, y = 0.55, hjust = 0, fontface = 2, label = pValue) 



pdf(file='ATA.pdf',width =3.5, height =3.5, onefile=F)
print(surPlot)
dev.off()



fit <- survfit(Surv(DSS, Event) ~Risk , data =HT0)  
diff=survdiff(Surv(DSS, Event) ~ Risk, data = HT0)  
pValue=1-pchisq(diff$chisq, df=1)
if(pValue<0.001){pValue=paste0("P= ",format(pValue, scientific = TRUE))}else{pValue=paste0("P= ",sprintf("%.03f",pValue))}    #pValue="p<0.001"
HR = (diff$obs[1]/diff$exp[1]) / (diff$obs[2]/diff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[1] + 1/diff$exp[2]))
HR = sprintf("%.02f",HR)
up95 = sprintf("%.02f",up95)
low95 = sprintf("%.02f",low95)
printHR = c(HR, low95, up95)
printHR 




pdf('HT0.pdf',width =3.5, height =3.5, onefile=F)
ggsurvplot(fit, 
           data = rt,   
           conf.int = F,  
           pval = pValue,
           pval.size= 5,
           pval.coord = c(0.2,0.7),
           ggtheme = theme_test(base_size = 12, base_line_size = 0.2,base_rect_size = 0.8) +
             theme(plot.title=element_text(hjust=0.5)),  
           title = paste0("HT0","HR = ",HR,"[",low95,"-",up95,"]"),   
           palette = mycol0,  
           legend = c(0.7, 0.3),  
           font.legend = 12,
           font.main = c(12, "plain","black"),
           font.x = c(12,"plain", "black"),
           font.y = c(12,"plain", "black"),
           font.tickslab = c(12, "plain", 'black'),
           xlab="Time(years)",
           ylab="TCFi probability",
           ylim = c(0.5,1),
           break.time.by = 2
) 




dev.off()  

