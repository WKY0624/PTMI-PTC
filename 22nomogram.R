library(foreign)
library(rms)
library(survival)

riskFile="totalRisk.txt"  
cliFile="06.clinical477.txt"   

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$RiskScore = as.numeric(risk$RiskScore)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('NA',x)))),,drop=F]

newcolnames=c("Tstage","Tstage1","Nstage","Mstage")
colnames(cli)[5:8]=newcolnames
cli$Tstage=as.factor(cli$Tstage)


samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
colnames(risk1)[1]="RFS"  
colnames(risk1)[2]="Recurrence"  
rt=cbind(risk1[,c("RFS", "Recurrence", "Risk")], cli)
rt=rt[,grep("RFS|Recurrence|Risk|Age|Gender|Tstage|Nstage|Mstage",colnames(rt))] 
rt<-na.omit(rt)
colnames(rt)[colnames(rt)=="Risk"] <- "RiskGroup"   
colnames(rt)[colnames(rt)=="Tstage"] <- "TStage"  
colnames(rt)[colnames(rt)=="Nstage"] <- "NStage"  
colnames(rt)[colnames(rt)=="Mstage"] <- "MStage" 

dd<-datadist(rt)
options(datadist='dd') 

coxm <-cph(Surv(RFS,Recurrence==1)~RiskGroup+Age+Gender+TStage+NStage+MStage, data=rt, x=T, y=T, surv=T)
surv<- Survival(coxm) 
surv1<- function(x)surv(1,lp=x) 
surv2<- function(x)surv(3,lp=x)
surv3<- function(x)surv(5,lp=x) 


plot(nomogram(coxm, fun=list(surv1,surv2,surv3),
              lp=F,  
              funlabel=c("1-Year RFS",'3-Year RFS','5-Year RFS'),
              maxscale=100,  
              fun.at=c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  ),
xfrac=.25, 
label.every = 1,  
col.grid = gray(c(0.8,0.95)),   

)   

dev.off()


nomoRisk=predict(coxm, data=rt, type="lp")
rt2=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt2), rt2)  
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)


f<-coxph(Surv(RFS,Recurrence==1)~Age+Gender+TStage+NStage+MStage+RiskGroup,data=rt)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index
write.table(c_index,file = "cindex.txt", sep = '\t', col.names = F,quote = F)


library(Hmisc)
S<-Surv(rt$RFS,rt$Recurrence==1)
rcorrcens(S~predict(coxm),outx=TRUE)


f <- cph(Surv(RFS, Recurrence) ~ Nomogram, x=T, y=T, surv=T, data=rt2, time.inc=1)
cal1 <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt2)/3), B=1000)

f <- cph(Surv(RFS, Recurrence) ~ Nomogram, x=T, y=T, surv=T, data=rt2, time.inc=3)
cal3 <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt2)/3), B=1000)
f <- cph(Surv(RFS, Recurrence) ~ Nomogram, x=T, y=T, surv=T, data=rt2, time.inc=5)
cal5 <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt2)/3), B=1000)

pdf(file="calibration.pdf", width=5, height=5)
plot(cal1, xlim=c(0.6,1), ylim=c(0.6,1),
     xlab="Nomogram-predicted RFS", 
     ylab="Observed RFS", 
     lwd=1,
     lty=5,
     col="green",
     sub=F)
plot(cal3, xlim=c(0.6,1), ylim=c(0.6,1), xlab="", ylab="", 
     lwd=1, lty=5, col="blue", sub=F, add=T)
plot(cal5, xlim=c(0.6,1), ylim=c(0.6,1), xlab="", ylab="",  
     lwd=1, lty=5, col="red", sub=F, add=T)

lines(cal5[,c("mean.predicted","KM")],type="b",lwd=3,col='red', pch=4)
lines(cal3[,c("mean.predicted","KM")],type="b",lwd=3,col='blue', pch=4)
lines(cal1[,c("mean.predicted","KM")],type="b",lwd=3,col='green', pch=4)
legend('bottomright', c('1-Year', '3-Year', '5-Year'),
       col=c("green","blue","red"), lwd=2, bty = 'n')

dev.off()




modelFile="model.txt"       

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$RiskScore = as.numeric(risk$RiskScore)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('NA',x)))),,drop=F]
newcolnames=c("Tstage","Nstage","Mstage")
colnames(cli)[6:8]=newcolnames
cli$Tstage=as.factor(cli$Tstage)  
model=read.table(modelFile,header=T, sep="\t", check.names=F, row.names=1)
model$Nomogram=as.numeric(model$Nomogram)
model$TNM=as.factor(model$TNM)   
model$ATA=as.factor(model$ATA)
model$MACIS=as.factor(model$MACIS)
model$EORTC=as.factor(model$EORTC) 


samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
model = model[samSample,,drop=F]
colnames(risk1)[1]="RFS"  
colnames(risk1)[2]="Recurrence" 
rt=cbind(risk1[,c("RFS", "Recurrence", "RiskScore")], model[,c("Nomogram")],cli)
colnames(rt)[4] <- "Nomogram"
rt<-na.omit(rt)

dd<-datadist(rt)
options(datadist='dd') 


f<-coxph(Surv(RFS,Recurrence==1)~ Nomogram, data=rt)   
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index  
cindex <- read.table("5Cindex.txt",header=T, sep="\t", check.names=F)



library(ggprism)
library(ggplot2)

pdf('Cindex.pdf', height = 3, width = 4)   
ggplot(cindex, aes(x=factor(model, levels = c(unique(cindex$model))), y=Cindex)) +
  geom_bar(stat = 'identity', 
           fill = 'white',   
           color = c('red','green2',"blue",'orange','purple'),   
           position= "dodge",
           width = 0.7  
  ) + 
  coord_cartesian(ylim = c(0.5, 0.8)) +
  geom_errorbar(aes(ymin=Cindex-SE, ymax=Cindex+SE), width=.2,
                color = c('red','green2',"blue",'orange','purple')) + 
  labs(y = 'C-index', x="") + 
  theme_prism(base_fontface = "plain",
              base_size = 12, 
              base_line_size = 0.5, 
              axis_text_angle = 45) +
  theme(plot.margin = margin(t = 10,  
                             r = 10,  
                             b = 0,  
                             l = 10))

dev.off()


library(survival)
library(dplyr)
rt_split <- survSplit(rt, cut = 1:5, end = "RFS", event = "Recurrence", start = "start", episode = "year")
c_index_vec <- numeric(5)
se_vec <- numeric(5)
for (t in 1:5) {
  f <- coxph(Surv(start, RFS, Recurrence) ~ MACIS, data = subset(rt_split, year == t))   
  sum_f <- summary(f)
  c_index_vec[t] <- sum_f$concordance[1]
  se_vec[t] <- sum_f$concordance[2]
}
result_mat <- rbind(c_index_vec, se_vec)
rownames(result_mat) <- c("C-index", "SE")
colnames(result_mat) <- paste0(1:5, " year")
result_mat <- as.data.frame(t(result_mat))
result_mat$model <- "Nomogram"    
result_mat$year <- c(1:5)
Nomogram <- result_mat       
result_mat$model <- "ATA"   
result_mat$year <- c(1:5)
ATA <- result_mat
result_mat$model <- "TNM" 
result_mat$year <- c(1:5)
TNM <- result_mat       
result_mat$model <- "MACIS"    
result_mat$year <- c(1:5)
MACIS <- result_mat       
result_mat$model <- "EORTC"    
result_mat$year <- c(1:5)
EORTC <- result_mat      

cindex <- rbind(Nomogram,ATA,TNM,MACIS,EORTC)
colnames(cindex)[1] <- "Cindex"


label_data <- cindex %>%
  group_by (model) %>%
  summarize(year = nth(year, 5), Cindex = nth(Cindex, 5), .groups = 'drop')

ggplot(cindex, aes(x = year, y = Cindex, color = model)) +
  geom_smooth(method = 'loess', formula = y ~ x, se = F) +
  geom_hline(yintercept = 0.5, lty = 2, color = "gray") +
  scale_x_continuous(limits=c(1,4), breaks=seq(1, 6, 1)) +        
  ylim(0.5,0.9) + 
  scale_color_manual(values = c('#e11a0c', '#48af45','#337cba','orange','purple'),
                     breaks = c('Nomogram','TNM',"ATA","MACIS","EORTC")) +
  labs(x = "Time (year)", y = "C-index") +
  theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme_classic(base_line_size = 0.3)+
  theme(axis.title.x=element_text(size = 10, vjust = -2),  
        axis.title.y=element_text(size = 10, vjust = 2),  
        axis.text.x=element_text(size = 10, vjust = -1), 
        axis.text.y=element_text(size = 10, hjust = -1),
        panel.grid.minor = element_blank())  + 
  theme(legend.text=element_text(size = 8),legend.spacing.y = unit(0,"cm"), legend.position = "none" ) +  #c(0.8,0.8)
  theme(plot.margin = margin(t = 15, r = 10, b = 30, l = 30)) + 
  geom_text(data=label_data, aes(x=year,y=Cindex,label=model), size=3.5, 
            hjust = c(-0.2,-0.1,-0.1, 0.2, -0.1),  
            vjust=c(-0.3, 0.6, 1, -2, 0.1))  

library(rms)
library(ggDCA)
library(survival)  
library(ggprism)

rm(list = ls()) 

riskFile = "totalRisk.txt"     
cliFile="06.clinical477.txt"  
modelFile="model.txt"     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$RiskScore = as.numeric(risk$RiskScore)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('NA',x)))),,drop=F]
newcolnames=c("Tstage","Nstage","Mstage")
colnames(cli)[6:8]=newcolnames
cli$Tstage=as.factor(cli$Tcategory2)   
model=read.table(modelFile,header=T, sep="\t", check.names=F, row.names=1)
model$Nomogram=as.numeric(model$Nomogram)
model$TNM=as.factor(model$TNM)  
model$ATA=as.factor(model$ATA)
model$MACIS=as.factor(model$MACIS)
model$EORTC=as.factor(model$EORTC) 

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
model=model[samSample,,drop=F]
colnames(risk1)[1]="RFS"  
colnames(risk1)[2]="Recurrence"  
rt=cbind(risk1[,c("RFS", "Recurrence", "RiskScore")], cli, model)
rt<-na.omit(rt)

dd<-datadist(rt) 
options(datadist='dd') 
Free <- coxph(Surv(RFS,Recurrence==1) ~ Age+Gender+Tstage+Nstage+Mstage+RiskScore, data=rt)
Score <- coxph(Surv(RFS,Recurrence==1) ~ RiskScore, data=rt)
Nomogram <- coxph(Surv(RFS,Recurrence==1) ~ Age+Gender+Tstage+Nstage+Mstage+RiskScore, data=rt)
ATA <- coxph(Surv(RFS,Recurrence==1) ~ ATA, data=rt)
TNM <- coxph(Surv(RFS,Recurrence==1) ~ TNM, data=rt)
MACIS <- coxph(Surv(RFS,Recurrence==1) ~ MACIS, data=rt)
EORTC <- coxph(Surv(RFS,Recurrence==1) ~ EORTC, data=rt)

P.ATA <- coxph(Surv(RFS,Recurrence==1) ~ RiskScore+ATA, data=rt)
P.TNM <- coxph(Surv(RFS,Recurrence==1) ~ RiskScore+TNM, data=rt)
P.MACIS <- coxph(Surv(RFS,Recurrence==1) ~ RiskScore+MACIS, data=rt)
P.EORTC <- coxph(Surv(RFS,Recurrence==1) ~ RiskScore+EORTC, data=rt)
dca_1<-dca(Nomogram,ATA,TNM,MACIS,EORTC, times=1)
dca_3<-dca(Nomogram,ATA,TNM,MACIS,EORTC, times=3)
dca_5<-dca(Nomogram,ATA,TNM,MACIS,EORTC, times=5)
dca_10 <- dca(Nomogram,ATA,TNM,MACIS,EORTC, times=10)
dca_13510<- dca(Nomogram,ATA,TNM,MACIS,EORTC, times=c(1,3,5,10))
dca_3510<- dca(Nomogram,ATA,TNM,MACIS,EORTC, times=c(3,5,10))
dca_35<- dca(Nomogram,ATA,TNM,MACIS,EORTC, times=c(3,5))
dca_135<- dca(Nomogram,ATA,TNM,MACIS,EORTC, times=c(1,3,5))

dca_Nomo1 <- dca(Nomogram, times=1)
dca_Nomo3 <- dca(Nomogram, times=3)
dca_Nomo5 <- dca(Nomogram, times=5)

dca_ATA1 <- dca(ATA, times=1)
dca_ATA3 <- dca(ATA, times=3)
dca_ATA5 <- dca(ATA, times=5)

pdf('DCA-5year.pdf', height = 4, width = 5.5)  

mycol=c('#e11a0c', '#337cba','#48af45','orange','purple','gray','black')

ggplot(dca_5,aes(x=thresholds, y=NB, group=model),
       lwd = 0.6) + 
  theme_classic(base_line_size = 0.3) +
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 1) +
  scale_x_continuous(limits = c(0,0.75), guide = 'prism_minor',breaks = seq(0, 1, 0.2)) +   
  scale_y_continuous(limits = c(-0.05,0.1), guide = 'prism_minor') +   
  scale_color_manual(values = c('red','blue','green','orange','purple','gray','black')) + 
   scale_linetype_manual(values = c('solid','solid','solid','solid','solid','twodash','twodash')) +
  theme(legend.title = element_blank() ,   
  legend.position = c(0.8,0.75),
  legend.text=element_text(size = 10),
  legend.spacing.y = unit(0,"cm"),
  legend.key.height = unit(0.5,"cm"))+
  annotate('text', x = 0.22, y = 0.03, label = "5-Year", colour="black", size = 5) +  
  theme_classic() +
  theme(axis.title.x=element_text(vjust=-5, size=15),  
        axis.title.y=element_text(size=15, vjust=5),  
        axis.text.x=element_text(size=5,
                                 vjust= -1),  
        axis.text.y=element_text(size=5, hjust = 0)) +  
  theme(plot.margin = margin(t = 15, r = 10, b = 30, l = 30))  

dev.off() 


library(survival)
library(timeROC)

riskFile="totalRisk.txt"  
cliFile="model.txt"        

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("DSS", "Event", "RiskScore")]   
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

rt$TNM = ifelse(rt$TNM=="Stage I",1,
                ifelse(rt$TNM=="Stage II",2,
                       ifelse(rt$TNM=="Stage III~IV",3,NA)))


rt$ATA = ifelse(rt$ATA=="Low",1,
                ifelse(rt$ATA=="Intermediate",2,
                       ifelse(rt$ATA=="High",3,NA)))

rt$MACIS = ifelse(rt$MACIS=="Group I~II",1,
                  ifelse(rt$MACIS=="Group III",2,
                         ifelse(rt$MACIS=="Group IV",3,NA)))
rt$EORTC = ifelse(rt$EORTC=="Group I~II",1,
                  ifelse(rt$EORTC=="Group III",2,
                         ifelse(rt$EORTC=="Group IV~V",3,NA)))


ROC_rt=timeROC(T=rt$DSS, delta=rt$Event,
               marker=rt$Nomogram, cause=1,    
               weighting='aalen',
               times=c(1,2,3,4,5,6,7,8,9,10), ROC=TRUE)

pdf(file="timeROC_Nomo.pdf", width=4.5, height=5)
par(tck = -0.02, 
    mgp = c(2, 0.5, 0))   
plot(ROC_rt,time=1,col='white',title=FALSE,lwd=2,
     cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=0.6, 
)
abline(h = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), v = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), col = gray(0.95,0.95))

plot(ROC_rt,time = 1,col=bioCol[1],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time = 3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time = 5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time = 10,col=bioCol[4],add=TRUE,title=FALSE,lwd=2)

legend('bottomright',   #x=0.4,y=0.3,
       c(paste0('  1-year (',sprintf("%.03f",ROC_rt$AUC[1]), ')'),
         paste0('  3-year (',sprintf("%.03f",ROC_rt$AUC[3]), ')'),
         paste0('  5-year (',sprintf("%.03f",ROC_rt$AUC[5]), ')'),
         paste0('10-year (',sprintf("%.03f",ROC_rt$AUC[10]), ')')),
      
       seg.len = 1, 
       x.intersp = 0.5,  
       col=bioCol[1:4], lty=1, lwd=2, bty = 'n') 
legend('topleft', 'AUC',bty='n')
dev.off()




predictTime=5  
aucText=c()

pdf(file="5modelROC_5year.pdf", width=4.5, height=5)

i=5
ROC_rt=timeROC(T=rt$DSS,
               delta=rt$Event,
               marker=rt$Nomogram, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)

par(tck = -0.02, 
    mgp = c(2, 0.5, 0))  
plot(ROC_rt, time=predictTime, col='white', title=FALSE, lwd=2,
     cex.lab=0.8, cex.axis=0.8, cex.main=1, cex.sub=0.6)
abline(h = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), v = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), col = gray(0.95,0.95))

for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$DSS,
                 delta=rt$Event,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol2[i-3], title=FALSE, lwd=2, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i]," (",sprintf("%.3f",ROC_rt$AUC[2]),")"))
}

legend("bottomright", aucText,
       col=bioCol2[1:(ncol(rt)-1)],
       lty = 1, lwd = 2, bty= "n", seg.len = 1, 
       x.intersp = 0.5,  
)

legend('topleft', '5-year AUC',bty='n')

dev.off()


ROC_rt_nomo=timeROC(T=rt$DSS, delta=rt$Event,
                    marker=rt$Nomogram, cause=1,
                    weighting='aalen',
                    times=c(1,2,3,4,5,6,7,8,9,10), ROC=TRUE)
Nomogram = as.data.frame(ROC_rt_nomo[["AUC"]])
colnames(Nomogram) = "AUC"
Nomogram$Group = "Nomogram"
Nomogram$Time = c(1:10)
Nomogram


ROC_rt_TNM=timeROC(T=rt$DSS, delta=rt$Event,
                   marker=rt$TNM, cause=1,
                   weighting='aalen',
                   times=c(1,2,3,4,5,6,7,8,9,10), ROC=TRUE)
TNM = as.data.frame(ROC_rt_TNM[["AUC"]])
colnames(TNM) = "AUC"
TNM$Group = "TNM"
TNM$Time = c(1:10)
TNM

ROC_rt_ATA=timeROC(T=rt$DSS, delta=rt$Event,
                   marker=rt$ATA, cause=1,
                   weighting='aalen',
                   times=c(1,2,3,4,5,6,7,8,9,10), ROC=TRUE)
ATA = as.data.frame(ROC_rt_ATA[["AUC"]])
colnames(ATA) = "AUC"
ATA$Group = "ATA"
ATA$Time = c(1:10)
ATA

ROC_rt_MACIS=timeROC(T=rt$DSS, delta=rt$Event,
                     marker=rt$MACIS, cause=1,
                     weighting='aalen',
                     times=c(1,2,3,4,5,6,7,8,9,10), ROC=TRUE)
MACIS = as.data.frame(ROC_rt_MACIS[["AUC"]])
colnames(MACIS) = "AUC"
MACIS$Group = "MACIS"
MACIS$Time = c(1:10)
MACIS

ROC_rt_EORTC=timeROC(T=rt$DSS, delta=rt$Event,
                     marker=rt$EORTC, cause=1,
                     weighting='aalen',
                     times=c(1,2,3,4,5,6,7,8,9,10), ROC=TRUE)
EORTC = as.data.frame(ROC_rt_EORTC[["AUC"]])
colnames(EORTC) = "AUC"
EORTC$Group = "EORTC"
EORTC$Time = c(1:10)
EORTC

AUCs = rbind(Nomogram,ATA,TNM,MACIS,EORTC)
write.table(AUCs, "AUCs.txt", row.names = F, col.names=T, sep="\t",quote=F)


library(ggplot2)
library(dplyr)

AUCdata=read.table('AUCs.txt', header=T, sep="\t", check.names=F)
label_data <- AUCdata %>%
  group_by(Group) %>%
  summarize(Time = Time[10], AUC = AUC[10],.groups="drop")

ggplot(AUCdata, aes(x = Time, y = AUC, color = Group)) +
  geom_smooth(method = 'loess', formula = y ~ x, se = F) +
  geom_hline(yintercept = 0.5, lty = 2, color = "gray") +
  
  scale_x_continuous(limits=c(1,12), breaks=seq(1, 10, 1)) +     
  ylim(0.55,0.85) + 
  scale_color_manual(values = c('#e11a0c', '#48af45','#337cba','orange','purple'),
                     breaks = c('Nomogram','TNM',"ATA","MACIS","EORTC")) +
  labs(x = "Time (year)", y = "AUC") +
  theme_classic(base_line_size = 0.3)+
  theme(axis.title.x=element_text(size = 10, vjust = -2), 
        axis.title.y=element_text(size = 10, vjust = 2),  
        axis.text.x=element_text(size = 10, vjust = -1), 
        axis.text.y=element_text(size = 10, hjust = -1), 
        panel.grid.minor = element_blank())  + 
  theme(legend.text=element_text(size = 8),legend.spacing.y = unit(0,"cm"), legend.position = "none" ) +  #c(0.8,0.8)
  theme(plot.margin = margin(t = 15, r = 10, b = 30, l = 30)) + 
  geom_text(data=label_data, aes(x=Time,y=AUC,label=Group), size=3.5, 
            hjust = c(-0.2,-0.1,-0.1, 0.2, -0.1),  
            vjust=c(-0.3, 0.6, 1, -2, 0.1))   

ggsave("AUCs.pdf", height = 3.5, width = 4)
dev.off()


