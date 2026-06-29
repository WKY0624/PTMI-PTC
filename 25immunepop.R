library(limma)
library(scales)
library(ggplot2)
library(ggtext)
library(reshape2)
library(tidyverse)
library(ggpubr)

riskFile="totalRisk.txt"    
immFile="infiltration_estimation_for_tcga.csv"   

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))  
immune=avereps(immune) 

sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "RiskScore"]
immune=immune[sameSample,]  

x=as.numeric(risk)
x[x>quantile(x,0.99)]=quantile(x,0.99)
outTab=data.frame()  
for(i in colnames(immune)){   
  y=as.numeric(immune[,i])
  if(sd(y)<0.001){next}   
  corT=cor.test(x, y, method="spearman")  
  cor=corT$estimate    
  pvalue=corT$p.value
  if(pvalue<0.05){   
    outTab=rbind(outTab,cbind(immune=i, cor, pvalue))    
  }}


write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)   
data <- outTab
data$Algorithm = sapply(strsplit(data[,1],"_"), '[', 2)
data$Algorithm = factor(data$Algorithm,level=as.character(unique(data$Algorithm[rev(order(as.character(data$Algorithm)))])))
data$immune = sapply(strsplit(data[,1],"_"),'[', 1)
ysort = unique(data$immune)
data = data[order(data$immune, decreasing = F),]   
colslabels = rep(hue_pal()(length(levels(data$Algorithm))),table(data$Algorithm)) 
data$cor = as.numeric(data$cor)
data$pvalue = as.numeric(data$pvalue)

ggplot(data, aes(x = immune, y = cor)) +
  geom_segment(aes(x = immune, xend = immune, y = 0, yend = cor),
               size = 0.3, linetype = 'solid' , color = 'gray30') +   
  geom_point(data = data, aes(size = -log10(pvalue), color = Algorithm), shape = 19, alpha =0.9,
       
  ) +      
  scale_size(range = c(2,5)) +
  geom_hline(yintercept = 0, linetype = 2, color = 'gray20', size = 0.3) +
  scale_x_discrete(limits=factor(sort(ysort, decreasing = T))) +
  labs(x = "Immune cell",y = 'Spearman correlation coefficient', title = '') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right", 
        legend.key.size = unit(8,'pt'),
        legend.title = element_text(size = 10, lineheight = 4),
        legend.text = element_text(size = 8, lineheight = 4),
        axis.ticks = element_line(linewidth = 0.3),
        axis.text.y = element_text(color = 'black',
                                   hjust = 1, 
                                   size = 8,
                                   lineheight = 2),
        axis.title = element_text(size = 10, color = 'black'),
  ) 

ggsave('ImmCellCorLollipop.pdf', height = 5, width = 5)

