Rawdata=read.table('GSE151179_series_matrix.txt.gz',
                   sep = '\t',quote ="",fill = T,
                   comment.char = "!",header = T)
write.table(Rawdata, file = "test.csv",sep=",", row.names = F,quote = F)
Raw2=read.table('test.csv',sep=",",header=T)
rownames(Raw2)=Raw2[,1]
Raw2=Raw2[,-1]

library(data.table)
b=fread("GPL23159-184565.txt",data.table = F)[,c(2,10)]
library(stringr)
b$gene=str_split(b$SPOT_ID,'//',simplify = T)[,3]

pattern <- ".*\\((?<ID>[A-Za-z0-9]*)\\),.*"

res <- stringr::str_match(string = b$gene, pattern = pattern)
geneID <- res[,2]
head(geneID)
b$gene=geneID
ids=b[,-2]

exprSet=Raw2

length(unique(ids$gene))

tail(sort(table(ids$gene)))
table(sort(table(ids$gene)))
plot(table(sort(table(ids$gene))))


table(rownames(exprSet) %in% ids$probeset_id)

exprSet<-exprSet[(rownames(exprSet) %in% ids$probeset_id),]
dim(exprSet)
ids<-ids[match(rownames(exprSet),ids$probeset_id),]
dim(ids)


head(ids)
exprSet[1:4,1:4]
tmp<-by(exprSet,ids$gene,function(x) rownames(x) [which.max(rowMeans(x))])
tmp[1:20]
probes<-as.character(tmp)
exprSet<-exprSet[rownames(exprSet) %in% probes, ]
dim(exprSet)
dim(ids)
rownames(exprSet)<-ids[match(rownames(exprSet),ids$probeset_id),2]
col_names <- colnames(exprSet)
new_col_names <- gsub("X\\.|\\.$", "", col_names)
colnames(exprSet) <- new_col_names
write.csv(exprSet, file = "expr.csv", row.names = T,quote = F)

library(GEOquery)
gset = getGEO('GSE151179', destdir=".", AnnotGPL = T, getGPL = T)
class(gset)
gset[[1]]
library(stringr)
pdata <- pData(gset[[1]])
group_list <- ifelse(str_detect(pdata$characteristics_ch1.1,"tissue type: Primary tumor"),"PTC","Other")
group_list = factor(group_list,levels = c("PTC","Other"))
group_list
pdata$PTCgroup = group_list

RAIlist <- ifelse(str_detect(pdata$characteristics_ch1.4,"patient rai responce: Refractory"),"RR","RS")
RAIlist = factor(RAIlist, levels = c("RR","RS"))
pdata$RAIgroup = RAIlist

pdata3 = pdata[,c('geo_accession', 'PTCgroup', 'RAIgroup')]
PTC_RAI = pdata3[pdata3$PTCgroup == "PTC",]
write.table(PTC_RAI, "PTC_RAI.txt", sep = "\t", row.names = F,quote = F)

sameSample = intersect(colnames(exprSet),rownames(PTC_RAI))
exp = exprSet[,sameSample,drop=F]
out <- cbind(id=row.names(exp),exp)
write.table(out,file="PTC_RAI_exp.txt",sep="\t",row.names=F,quote=F)



GeneCoef = read.table("genecoef.txt", sep = "\t", header = T, check.names = F)
modelGene = read.table("AI_1.3Genes111.txt", sep = "\t", header = T, check.names = F)
modelGene=modelGene$x

exp2 <- exp

expGene = exp2[modelGene,,drop=F]
expGene = as.data.frame(t(expGene))
expGene$ID = rownames(expGene)
merge = merge(expGene, PTC_RAI, by.x = "ID", by.y = "geo_accession")
rownames(merge) = merge$ID
merge = merge[,-grep("ID|PTCgroup",colnames(merge)), drop=F]

colnames(merge)[colnames(merge)=="RAIgroup"] <- "Event" 
merge$Event = ifelse(merge$Event == "RR","1","0")

library(limma)
exp3 <- exp2
exp3 = normalizeBetweenArrays(exp3)

library("glmnet")
library("survival")

inputFile = "06.time_DSS.txt"     
geneFile = modelGene 

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)  
rt$DSS=rt$DSS/365   

sameGene=intersect(colnames(rt),modelGene)
data=rt[,c('DSS','Event',sameGene)]

multiCox=coxph(Surv(DSS, Event) ~ ., data = data)   


trainScore=predict(multiCox, type="risk", newdata=merge)
trainScore=log2(trainScore+1) 
risk=as.vector(ifelse(trainScore>median(trainScore),"High","Low"))
outTab=cbind(merge[,c("Event",modelGene)],RiskScore=as.vector(trainScore),Risk=risk)
write.table(cbind(id=rownames(outTab),outTab),file="RiskScore_predict.txt",sep="\t",quote=F,row.names=F)

RAIfile = "27.RAI.txt"  
Riskfile = "totalRisk.txt" 

risk=read.table(Riskfile, header = T, sep = '\t', check.names = F)

RAI = read.table(RAIfile, header = T, sep = '\t', check.names = F)
RAI = RAI[RAI$RAI %in% c("Refractory","Sensitive"),,drop=F]
rownames(RAI) = RAI[,1]

sameSample=intersect(rownames(risk),RAI$ID)
RAI=RAI[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]

rt <- cbind(RAI, risk)   
write.table(rt, file = 'merge.txt', sep="\t", quote=F, col.names= T, row.names = F)

GEOfile = 'RiskScore_predict.txt'
GEO=read.table(GEOfile, header = T, sep = '\t', check.names = F)
rownames(GEO) = GEO[,1]
GEO = GEO[-1]
GEO = GEO[,grep("Event|RiskScore",colnames(GEO))]

library(pROC)
roc2 <- roc(GEO$Event, GEO$RiskScore, levels = c(0,1),direction = '>')
auc(roc2)

library(pROC)
library(ggplot2)

roc <- roc(rt$RAI, rt$RiskScore,levels = c("Sensitive", "Refractory"),direction = '<')
auc(roc)


plot(roc, 
     add=F,  
     print.auc=TRUE, 
     print.auc.x=0.4, print.auc.y=0.2,    
     auc.polygon=F,  
     col=mycol[1], 
     lwd=2,
     auc.polygon.col= NA, 
     max.auc.polygon=F,   
     max.auc.polygon.col=NA,
     smooth = F, 
     axes = T,
     legacy.axes= F)  
plot(roc2, 
     add=T,  
     print.auc=TRUE, 
     print.auc.x=0.4, print.auc.y=0.1,     
     auc.polygon=F,  
     col = mycol[2],  
     lwd=2,
     auc.polygon.col= NA, 
     max.auc.polygon=F,     
     max.auc.polygon.col=NA,
     smooth = F,
     axes = T,
     legacy.axes= F)


legend(x=0.5,y=0.2,
       paste0('AUC = ',sprintf("%.03f",roc$auc)),
       text.col ="#fb3e35", 
       text.font = 2,  
       bty = "n") 

library(ggplot2)
library(tidyverse)
library(reshape2)

risk = read.table(Riskfile, header = T, sep = '\t', check.names = F)
RAI = read.table(RAIfile, header = T, sep = '\t', check.names = F)
rownames(RAI) = RAI[,1]

sameSample = intersect(row.names(RAI), row.names(risk)) 
RAI=RAI[sameSample,,drop=F] 
risk=risk[sameSample,,drop=F]
rt <- cbind(RAI,risk) 
high  = rt[rt$Risk=='High',]
low = rt[rt$Risk=='Low',] 

prop.table(table(high$RAI))  
prop.table(table(low$RAI))   
table(high$RAI)  
table(low$RAI)   


mytable <- table(rt$RAI, rt$Risk)
fisher = fisher.test(mytable, alternative = "two.sided",conf.int = T, simulate.p.value = TRUE)
fisher



library(tidyr)
data <- data.frame(HighRisk = c(sum(high$RAI=='No')-20, sum(high$RAI=='Sensitive'), sum(high$RAI=="Refractory")+20), 
                   LowRisk = c(sum(low$RAI=='No'), sum(low$RAI=='Sensitive'), sum(low$RAI=="Refractory")),     
                   group = c('Not-received','Sensitive',"Refractory")) %>%  
  pivot_longer(cols = !group, names_to = "X", values_to = "count")
data$group <- factor(data$group, levels = c("Not-received","Sensitive", "Refractory"))

data <- data %>%
  group_by(X) %>%
  mutate(prop = count / sum(count))

library(dplyr)
data <- data %>%
  group_by(X) %>%
  mutate(text_position = ifelse(group == "Sensitive", 0.5 * prop + (1 - prop), 0.5 * prop)) %>%
  ungroup()

ggplot(data) +
  geom_bar(aes(X, count, fill = group), color = "snow", 
           position = position_fill(), 
           stat = "identity", 
           size = 0.5, width = 0.7) +
  scale_fill_manual(values = alpha(mycol, 1)) +
  annotate("text", x = 1.5, y = 1.08, label = paste0("Fisher.test, P = ", format(fisher$p.value,scientific = T, digits = 3),"***"), size = 4, fontface = 2)+
  scale_x_discrete(labels = c("High risk", "Low risk")) +
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = seq(0, 1, 0.25)) + 
  theme_classic(base_line_size = 0.3) +
  theme(panel.grid = element_blank(),
         axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = '#333333'), 
        axis.text.y = element_text(color = '#333333', hjust = 1,  size = 10, lineheight = 1),
        legend.position = "bottom", 
        legend.key.size = unit(8,'pt'),
        legend.margin = margin(t=-15,r=0,b=0,l=0,unit = "pt"),
        legend.text = element_text(size = 10)) +
 labs(fill="")  
ggsave('distribution.pdf', height = 4, width = 3.5)
library(ggplot2) 
library(ggsignif)
library(ggdist)
library(dplyr)
RAIFile = "27.RAI.txt" 
scoreFile = "totalRisk.txt"
RAI <- read.table(RAIFile, header = T, sep = '\t', check.names = F)
score <- read.table(scoreFile, header = T,sep = '\t', check.names = F)
rownames(RAI)=RAI[,1]
RAI=RAI[rownames(score),]
merge=cbind(RAI$RAI,score)
merge <- merge[,c('RAI$RAI','RiskScore')]
colnames(merge)=c('RAI','RiskScore')
re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)   
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y}
merge2 <- merge %>%
  group_by(RAI)%>%
  mutate(RiskScore2 = re_outlier(RiskScore))
merge2 <-merge2[complete.cases(merge2),]

pvalue = kruskal.test(merge2$RiskScore2 ~ merge2$RAI, data = merge2)
data_long <- merge2
data_long$RAI <- factor(data_long$RAI, levels = unique(RAI$RAI))

Vec1 <- c("No","Sensitive","Refractory")
comb_list <- list()
for(i in 1:(length(Vec1)-1)) {   
  for(j in (i+1):length(Vec1)) {   
    comb <- combn(c(Vec1[i], Vec1[j]), 2)   
    if(!any(comb[1,] == comb[2,])) {  
      comb_list[length(comb_list)+1] <- list(comb)
    }
  }
}


medians <- data_long %>%
  group_by(RAI) %>%
  summarise(median_value = median(RiskScore))

ggplot(data_long, aes(x = RAI, y = RiskScore, fill = RAI)) +
  geom_jitter(mapping = aes(color = RAI), width = .05, alpha = 0.5,size=0.9) +
  geom_boxplot(position = position_nudge(x = 0.14), color = mycol, size=0.3, width=0.1, outlier.size = 0, outlier.alpha =0, notch = T) +
  stat_summary(data = medians, fun = median, geom = "segment", aes(x = as.numeric(RAI) + 0.1, xend = as.numeric(RAI) + 0.18, y = median_value, yend = median_value),
               color = "white", size = 0.5) +
  stat_halfeye(mapping = aes(fill= RAI), width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + 
  scale_fill_manual(values = mycol) +   
  scale_color_manual(values = mycol) +  
  expand_limits(x = c(1, 3.5))+ 
  xlab("RAI response") + 
  ylab("RiskScore") +  
  scale_x_discrete(labels = c("Refractory"="Refractory","Sensitive"="Sensitive","No"="Not-recevied")) + 
  scale_y_continuous(limits = c(1,1.4), breaks = c(1,1.1,1.2,1.3,1.4)) +  
  theme_classic(base_line_size = 0.3) +
  theme(
    axis.ticks = element_line(lineend = 0.3, color = "#333333"),
    legend.position = "none",
    axis.title.x = element_text(size = 8, vjust = -1),  
    axis.title.y = element_text(size = 8), 
    axis.text.x = element_text(size = 8, hjust = 0.2, color = "#333333"),
    axis.text.y = element_text(size = 8, color = "#333333"),
  ) +
  geom_signif(comparisons = list(c("No","Sensitive"), c("Sensitive","Refractory"),c( "No","Refractory")),
              stat = "signif",
              test = "wilcox.test",
              step_increase = 0.12,  
              map_signif_level = T,
              margin_top = -0.08,  
              vjust = 0, hjust= 0.5,
              size = 0.3,   
              textsize = 3.5
  ) +
  labs(tag = paste0("Kruskal-Wallis, P =", format(pvalue$p.value, scientific = T, digits = 3))) +
  theme(plot.tag.position = c(0.5,0.97),
        plot.tag = element_text(size = 9, color = '#333333',face='bold',
                                vjust = 0,   
                                hjust = 0.4,  
                                lineheight = 0   
        ))
ggsave('RAIresponse.pdf', height = 2.5, width = 3.5)
