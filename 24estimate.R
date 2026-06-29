library(limma)
library(estimate)
inputFile="05.TPM100.txt"    

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=log2(data+1)   
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]

out=rbind(ID=colnames(data),data)
write.table(out,file="uniq.symbol_TPM.txt",sep="\t",quote=F,col.names=F)

filterCommonGenes(input.f="uniq.symbol_TPM.txt", 
                  output.f="commonGenes_TPM.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes_TPM.gct",
              output.ds="estimateScore_TPM.gct")

scores=read.table("estimateScore_TPM.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)

write.table(out, file="TMEscores_TPM.txt", sep="\t", quote=F, col.names=F)

library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)
library(rstatix)  
library(gghalves)   

riskFile="totalRisk.txt" 
scoreFile="TMEscores_TPM.txt"    

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$Risk=factor(risk$Risk, levels=c("High","Low"))

score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score=as.matrix(score)
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))  
score=avereps(score)
score=score[,1:4]  
risk2=risk[row.names(score),"Risk",drop=F]   
score=score[row.names(score),,drop=F]  
rt=cbind(risk2, score)  
rt=rt[,c(1:5)]   


re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)   
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y}

df_is <- rt%>%
  group_by(Risk)%>%
  mutate(ImmuneScore = re_outlier(ImmuneScore))
df_is<-df_is[complete.cases(df_is),]

sam <- intersect(rt$TumorPurity, df_is$TumorPurity)
newdf<-rt[which(rt$TumorPurity%in%sam),]
newrow <- rownames(newdf)

data_all=melt(df_is, id.vars=c("Risk"))   
colnames(data_all)=c("RiskGroup", "scoreType", "Score")


data <- data_all[data_all$scoreType == 'StromalScore'|
                   data_all$scoreType == 'ImmuneScore'| 
                   data_all$scoreType == 'ESTIMATEScore',]
data$scoreType = ifelse(data$scoreType=='StromalScore','Stromal score',
                        ifelse(data$scoreType=='ImmuneScore','Immune score','ESTIMATE score'))  
data_tumor <- data_all[data_all$scoreType =='TumorPurity',]
data_tumor$scoreType = ifelse(data_tumor$scoreType=='TumorPurity','Tumor purity')
data_4in1 <- data_all[data_all$scoreType == 'StromalScore'|
                        data_all$scoreType == 'ImmuneScore'| 
                        data_all$scoreType == 'ESTIMATEScore'|
                        data_all$scoreType == 'TumorPurity',]
data_4in1$scoreType = ifelse(data_4in1$scoreType=='StromalScore','Stromal score',
                             ifelse(data_4in1$scoreType=='ImmuneScore','Immune score',
                                    ifelse(data_4in1$scoreType=='ESTIMATEScore','ESTIMATE score','Tumor purity')))   ##修改单元格名称————(4合1)
colnames(data_4in1)[colnames(data_4in1)=="RiskGroup"] <- "Risk group"   
data <- data %>%  group_by(scoreType)
data$Score <- as.numeric(data$Score)
stat.test_data <- wilcox_test(data = data,   
                              Score ~ RiskGroup, 
                              paired = F,   
                              p.adjust.method = "fdr",
) %>%
  add_xy_position(x = "scoreType")
stat.test_data$p.scient <- format(stat.test_data$p, scientific = TRUE)
stat.test_data$p.round3 <- round(stat.test_data$p, 3)
stat.test_data$p.adj.signif <- ifelse(stat.test_data$p<0.001,"***",ifelse(stat.test_data$p<0.01,"**",ifelse(stat.test_data$p<0.05,"*","")))
stat.test_data$p.value <- ifelse(stat.test_data$p<0.001, format(stat.test_data$p, scientific = TRUE),round(stat.test_data$p, 3))
data_tumor <- data_tumor %>%  group_by(scoreType)
data_tumor$Score <- as.numeric(data_tumor$Score)
stat.test_data_tumor <- wilcox_test(data = data_tumor,
                                    Score ~ RiskGroup,
                                    paired = F,   
                                    p.adjust.method = "fdr",
) %>% add_xy_position(x = "scoreType")
stat.test_data_tumor$p.scient <- format(stat.test_data_tumor$p, scientific = TRUE)
stat.test_data_tumor$p.round3 <- round(stat.test_data_tumor$p,3)
stat.test_data_tumor$p.adj.signif <- ifelse(stat.test_data_tumor$p<0.001,"***",ifelse(stat.test_data_tumor$p<0.01,"**",ifelse(stat.test_data_tumor$p<0.05,"*","")))
stat.test_data_tumor$p.value <- ifelse(stat.test_data_tumor$p<0.001, format(stat.test_data_tumor$p, scientific = TRUE),round(stat.test_data_tumor$p, 3))


p1 <- ggplot()+
  geom_half_violin(data = data %>% filter(RiskGroup == "High"),
                   aes(x = scoreType, y = Score), side = "l",size= 1,
                   colour="white", fill=mycol[1], alpha = 0.2, width = 1,
                   position = position_dodge(width = 0.2)
  ) +
  geom_half_violin(data = data %>% filter(RiskGroup == "Low"),
                   aes(x = scoreType,y = Score), side = "r", size= 1,
                   colour="white", fill=mycol[2], alpha = 0.2, width = 1,
                   position = position_dodge(width = 0.2)
  ) +
  geom_half_boxplot(data = data %>% filter(RiskGroup == "High"),
                    aes(x = scoreType, y = Score), width = 0.2, 
                    colour=mycol[1], outlier.shape = 4, outlier.size = 1,
                    fill=mycol[1],side = "l", alpha = 0.5, nudge = 0.01,errorbar.draw = F,
                    position = position_dodge(width = 1))+
  geom_half_boxplot(data = data %>% filter(RiskGroup == "Low"),
                    aes(x = scoreType, y = Score), width = 0.2, 
                    colour=mycol[2], outlier.shape = 4, outlier.size = 1,
                    fill=mycol[2],side = "r", alpha = 0.6, nudge = 0.01, errorbar.draw = F,
                    position = position_dodge(width = 1)
  ) +
  
  geom_line(data = data, aes(x = scoreType, y = Score, color = RiskGroup),
            stat = 'summary', fun=median, lty =1, size =2,
            position = position_dodge(width = 0.1)) +
  scale_color_manual(values = mycol,name="Group",labels = c("High-risk","Low-risk"))+
  stat_pvalue_manual(
    stat.test_data, 
    label = 'P = {p.value}\n{p.adj.signif}', 
    bracket.size = 0.3, 
    bracket.shorten = 0.15,  
    tip.length = 0.01,  
    size = 3.5
  ) + 
  xlab("") +
  ylab("Tumor microenvironment score") +
  scale_y_continuous(limits = c(-4000,4300)) +   
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5)+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black', margin = margin(0.2,0,0,0, 'cm')),    #x轴标签
    axis.ticks =element_line(linewidth = 0.3),
    axis.text.y = element_text(color = 'gray30', hjust = 1, # 左对齐
                               size = 9, lineheight = 1),
    plot.title = element_text(hjust = 0.5, size = 10),
    text = element_text(family = ""),
    legend.position = "top", 
    legend.key.size = unit(10,'pt'), 
    legend.justification = "centre")


p2 <- ggplot()+
  geom_half_violin(data = data_tumor %>% filter(RiskGroup == "High"),
                   aes(x = scoreType, y = Score), side = "l",linewidth= 1,
                   colour="white", fill=mycol[1], alpha = 0.2, width = 1,
                   position = position_dodge(width = 0.2)) +
  geom_half_violin(data = data_tumor %>% filter(RiskGroup == "Low"),
                   aes(x = scoreType,y = Score), side = "r", linewidth= 1,
                   colour="white", fill=mycol[2], alpha = 0.2, width = 1,
                   position = position_dodge(width = 0.2)) +
  geom_half_boxplot(data = data_tumor %>% filter(RiskGroup == "High"),
                    aes(x = scoreType, y = Score), width = 0.2, 
                    colour=mycol[1], outlier.shape = 4, outlier.size = 1,
                    fill=mycol[1],side = "l", alpha = 0.5, nudge = 0.01,errorbar.draw = F,
                    position = position_dodge(width = 1))+
  geom_half_boxplot(data = data_tumor %>% filter(RiskGroup == "Low"),
                    aes(x = scoreType, y = Score), width = 0.2, 
                    colour=mycol[2], outlier.shape = 4, outlier.size = 1,
                    fill=mycol[2],side = "r", alpha = 0.6, nudge = 0.01,errorbar.draw = F,
                    position = position_dodge(width = 1)) +

  geom_line(data = data_tumor, aes(x = scoreType, y = Score, color = RiskGroup),
            stat = 'summary', fun=median, lty =1, linewidth =2,
            position = position_dodge(width = 0.1)) +
  scale_color_manual(values = mycol,name="Group",labels = c("High-risk","Low-risk"))+
  stat_pvalue_manual(
    stat.test_data_tumor, 
    label = 'P = {p.value}\n{p.adj.signif}', 
    bracket.size = 0.3, 
    bracket.shorten = 0.15,  
    tip.length = 0.01,  
    y.position = 1,
    linewidth = 3.5
  ) + 
  xlab("") +
  ylab("Percentage") +
  scale_y_continuous(position = "right", limits = c(0.45,1.1)) +   
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5)+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black', margin = margin(0.2,0,0,0, 'cm')),    #x轴标签
    axis.ticks =element_line(linewidth = 0.3),
    axis.text.y = element_text(color = "gray30", size = 9, lineheight = 1, hjust = 1), # 左对齐
    axis.title.y.right = element_text(vjust = 2),
    plot.title = element_text(hjust = 0.5, size = 10),
    text = element_text(family = ""),
    legend.position = 'top',  
    legend.key.size = unit(10,'pt'),
    legend.justification = "centre") 

library(patchwork)
p1 + p2 + 
  plot_annotation(#title = "Wilcox.test,
  ) +
  plot_layout(ncol=2,   
              widths = c(3, 1), 
              guides='collect'  
  ) & theme(legend.position='top',
            legend.key.size = unit(0.3,"cm"),
            legend.text = element_text(size = 10)) 

ggsave("4score.pdf", height = 3.5, width = 6) 



library(ggplot2)
library(gghalves)
library(dplyr)
library(ggsignif)
library(ggpubr)

dt=rt
dt <- as_tibble(dt)
head(dt)


p <- ggplot(dt, aes(x = Cluster, y = TumorPurity, 
                    color = Cluster, fill = Cluster))+
  stat_summary(fun = mean, 
               geom = 'bar', 
               width = 0.4,
               alpha = 0.8,
               show.legend = F) +
  scale_fill_manual(values = mycol) + 
  scale_color_manual(values = mycol) + 
  theme_minimal()
p
p1 <- p +
  stat_summary(fun = mean,
               geom = "errorbar",
               width = 0.3,
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),
               show.legend = F)
p1
p2 <- p1 +
  geom_half_violin(side = 'r', 
                   nudge = 0.3, 
                   scale = 'width',
                   width = 0.5,
                   adjust = 0.8,
                   trim = F,
                   color = NA,
                   alpha = 0.6,
                   show.legend = F)
p2
p3 <- p2 +
  geom_half_point_panel(side = 'r', 
                        shape = 21, 
                        size = 1, 
                        color = 'white',
                        show.legend = F)
p3
p4 <- p3 + coord_flip() 
p4

p5 <- p4 +
  geom_signif(comparisons = list(c('C1', 'C2'), 
                                 c('C1', 'C3'),
                                 c('C2', 'C3')), 
              test = 'wilcox.test',
              y_position = c(max(dt$TumorPurity),
                             max(dt$TumorPurity) + 0.2,
                             max(dt$TumorPurity) + 0.4), 
              map_signif_level = T,
              color = 'black',
              show.legend = F)
p5

ggsave('tumor.pdf',height = 4,width = 6)


riskScore = risk[rownames(score),"RiskScore",drop=F]
riskScore = cbind(riskScore,score)

re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)  
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y}

riskScore2 <- riskScore%>%
  mutate(RiskScore = re_outlier(RiskScore))
riskScore2 <-riskScore2[complete.cases(riskScore2),]

SpearmanR <- cor(riskScore2$RiskScore, riskScore2$StromalScore,method="spearman",use="complete.obs")
SpearmanP <- cor.test(riskScore2$RiskScore, riskScore2$StromalScore, method="spearman", use="complete.obs")
ggplot(riskScore2, aes(x=RiskScore, y=StromalScore)) +
  geom_point(color = mycol2[4], alpha=0.7, pch=20, size=2) +
  geom_smooth(method=lm , formula = y ~ x, 
            
              color=mycol2[4], fill=mycol2[4], alpha = 0.3, se=TRUE) +
  theme_test(base_line_size = 0.3)+
  ylab("StromalScore") + 
  xlab("Score") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 10),
    axis.title.y = element_text(vjust = 0),
    axis.text = element_text(color = "gray30",size = 9)
  ) +   
  xlim(min(riskScore2$RiskScore), max(riskScore2$RiskScore)) + 
  ylim(min(riskScore2$StromalScore), max(riskScore2$StromalScore))+  
  annotate("text", x = max(riskScore2$RiskScore), y = max(riskScore2$StromalScore-400), fontface = 1, hjust = 1, label = paste0("rho = ",round(SpearmanR,3), "\n","P = ", format(SpearmanP$p.value,scientific = TRUE)))
ggsave("StromalScore.pdf", height = 3, width = 3.2)

