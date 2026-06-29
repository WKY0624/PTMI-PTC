library(survival)   

riskFile = "totalRisk.txt"
cliFile = "06.clinical477.txt" 

risk = read.table(riskFile, header=T, sep="\t", check.names=F,row.names=1)    
cli = read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)  
newcolnames=c("Tstage","T","Nstage","Mstage")
colnames(cli)[5:8]=newcolnames

cli = cli[,grep("Age|Gender|RadiationHistory|Tstage|Nstage|Mstage|Multifocality|BRAF|RAS|TERT|Genefusion|RiskScore",colnames(cli))]  
cli$Tstage=as.factor(cli$Tstage)

risk$RiskScore=scale(risk$RiskScore)

sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(RFS=risk[,1],    
         Recurrence=risk[,2], cli, RiskScore=risk[,(ncol(risk)-1)])
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(RFS, Recurrence) ~ rt[,i], data = rt)  
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
uniTab
write.table(uniTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

rt1=rt[,c("RFS", "Recurrence", as.vector(uniTab[,"id"]))]   
multiCox = coxph(Surv(RFS, Recurrence) ~ ., data = rt1)   
multiCoxSum = summary(multiCox)
multiTab = data.frame()
multiTab = cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=uniTab$id,multiTab)
multiTab
write.table(multiTab,file="multiCoxP.txt",sep="\t",row.names=F,quote=F)


library(forplo)

uni = read.table('uniCox.txt', header = T, row.names = 1, check.names = F, sep = '\t') 
uni2 = uni[,1:3]  
pdf('uniCox.pdf', height = 8, width = 15)
forplo(uni2,    
       em = "HR",
       row.labels = c('Age','Gender','RadiationHistory', 'T Stage','N Stage','M Stage','Multifocality','BRAF','RAS','TERT','RiskScore'),
       pval = sprintf("%.03f",uni$pvalue),    
       xlim = c(0.05,20),
       ci.sep = '-',
       ci.lwd = 2,
       ci.edge = T,   
       char = 17,   
       size = 1.2,   
       insig.col = '#d5dedd',
       fill.colors = c("#5b679b","#d5dedd","#d5dedd","#5b679b","#d5dedd","#5b679b","#d5dedd","#d5dedd","#d5dedd","#5b679b","#5b679b"),
       right.bar = T,   
       rightbar.ticks = T,  
       left.bar = T,
       leftbar.ticks = T,
       left.align = F,  
       margin.top = 10,
       margin.bottom = 5,
       margin.right = 1,
       title = ' ')
dev.off()

multi = read.table('multiCoxP.txt', header = T, row.names = 1, check.names = F, sep = '\t')
multi2 = multi[,1:3]
pdf('multiCoxP.pdf', height = 8, width = 15)
forplo(multi2,
       em = "HR",
       pval = sprintf("%.03f",multi$pvalue),
       xlim = c(0.05,50),
       ci.sep = '-',
       ci.lwd = 2,
       ci.edge = T,   
       char = 16,  
       size = 1.2,  
       insig.col = '#d5dedd',
       right.bar = T,   
       rightbar.ticks = T,  
       left.bar = T,
       leftbar.ticks = T,
       left.align = F,  
       margin.top = 10,
       margin.bottom = 5,
       margin.right = 1,
       title = ' ')
dev.off()

