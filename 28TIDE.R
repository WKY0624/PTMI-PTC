library(limma)
library(sva)
tcgaExpFile ="total.normalize.txt"

rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]   
exp=rt[,2:ncol(rt)]   
dimnames=list(rownames(exp),colnames(exp))  
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga2 <- t(apply(tcga, 1, function(x)x-(mean(x))))

tcgaTab=rbind(ID=colnames(tcga2), tcga2)
write.table(tcgaTab, file="TIDE_total.normalize", sep="\t", quote=F, col.names=F)

TIDEfile="TIDE.csv"  
Riskfile = "totalRisk.txt"  

tide=read.table(TIDEfile, header=T, sep=",", check.names=F) 
risk=read.table(Riskfile, header = T, sep = '\t', check.names = F)
rownames(tide) = tide[,1]
sameSample = intersect(row.names(tide), row.names(risk))
tide = tide[sameSample,,drop=F]
risk = risk[sameSample,,drop=F]
rt =cbind(tide,risk)  
write.table(rt, file = 'merge.txt', sep="\t", quote=F, col.names= T, row.names = F)

library(pROC)
library(ggplot2)

roc <- roc(rt$Responder, rt$RiskScore)
auc(roc)

plot(roc, 
     print.auc=TRUE, 
     print.auc.x=0.4, print.auc.y=0.2,   
     auc.polygon=T,  
     col="#fb3e35", 
     auc.polygon.col="#fff7f7",   
     max.auc.polygon=FALSE,   
     smooth = F,   
     legacy.axes= F)  


library(magrittr)
library(tidyr)
library(dplyr)

rt = cbind(tide, risk[,c(15:16)])  
high  = rt[rt$Risk=='High',]
low = rt[rt$Risk=='Low',] 

prop.table(table(high$Responder))  
prop.table(table(low$Responder))   
table(high$Responder)  
table(low$Responder)   

rt$Responder <- as.logical(rt$Responder)
ks.test(scale(rt$Responder),'pnorm')  
mytable <- table(rt$Responder, rt$Risk)

chisq.test(mytable)
chisq=chisq.test(mytable)

data <- data.frame(HighRisk = c(sum(high$Responder=='False'), sum(high$Responder=='True')), 
                   LowRisk = c(sum(low$Responder=='False'), sum(low$Responder=='True')),    
                   group = c('Non-Responder','Responder')) %>%  
  pivot_longer(cols = !group, names_to = "X", values_to = "count")
data$group = sort(data$group, decreasing = F)

data <- data %>%
  group_by(X) %>%
  mutate(prop = count / sum(count))

data <- data %>%
  group_by(X) %>%
  mutate(text_position = ifelse(group == "Non-Responder", 0.5 * prop + (1 - prop), 0.5 * prop)) %>%
  ungroup()


ggplot(data) +
  geom_bar(aes(X, count, fill = group), color = "#f3f4f4",
           position = position_fill(), 
           stat = "identity", 
           size = 0.5, width = 0.6) +
  scale_fill_manual(values = alpha(mycol, 0.9), labels = rev(c('Respond', 'Non-Respond'))) +
  geom_text(aes(X, text_position, label = paste0(round(prop * 100), "%")), 
            size = 4, color = "white", fontface = 2) +
  annotate("text", x = 1.5, y = 1.08, label = paste0("chisq.test, P = ", round(chisq$p.value,3)), size = 4, fontface = 2)+
  scale_x_discrete(labels = c("High risk", "Low risk")) +
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = seq(0, 1, 0.25)) + 
  
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = '#333333'), 
        axis.text.y = element_text(color = '#333333', hjust = 1,  size = 10, lineheight = 1),
        legend.position = "top",
        legend.key.size = unit(10,'pt'),
        legend.text = element_text(size = 10)) +
  guides(fill=guide_legend(reverse=T)) + 
  labs(fill="")  

ggsave('distribution.pdf', height = 4, width = 3)


library(limma)
library(reshape2)
library(ggpubr)
library(ggridges)

library(rstatix)  
library(ggpubr)

data_multi <- rt[,c("Patient",'Risk','TIDE','IFNG','MSI Expr Sig','Merck18','CD274','CD8','Dysfunction','Exclusion','MDSC','CAF','TAM M2')]
rownames(data_multi) = data_multi$Patient
data_multi = data_multi[,-1]
data_multi$MDSC = data_multi$MDSC*10
data_multi$CAF = data_multi$CAF*10
data_multi$`TAM M2` = data_multi$`TAM M2`*10


ks.test(scale(rt$TIDE),'pnorm')   
ks.test(scale(low$TIDE),'pnorm')  


wilcox.test(CD274 ~ Risk, data = data_multi, alternative = 'two.sided')

stat_multi <- data.frame(Tides=colnames(data_multi)[2:ncol(data_multi)])
for (i in 2:ncol(data_multi)){
  print(i)
  stat_multi[i-1,2] <- wilcox.test(data_multi[,i] ~ Risk, data = data_multi, 
                                   alternative = 'two.sided',
                                   exact = FALSE)[["p.value"]]}
stat_multi$fdr <- p.adjust(stat_multi$V2, method = "fdr")
stat_multi$bon <- p.adjust(stat_multi$V2, method = 'bonferroni')

colnames(stat_multi) <- c("Tides","P.value","FDR","Bonferroni")

stat_multi$P.round3 <- round(stat_multi$P.value, 3)
stat_multi$P.sig <- ifelse(stat_multi$P.value<0.001,"***",ifelse(stat_multi$P.value<0.01,"**",ifelse(stat_multi$P.value<0.05,"*","")))


data2 = reshape2::melt(data_multi, id.vars = "Risk")
colnames(data2) = c("Risk","Pathway","Score")

ggboxplot(data2, x='Pathway', y= "Score", color = "Risk",
          notch = F, size = 0.4, width = 0.6, outlier.shape = 3, outlier.size = 1,
          xlab="",ylab="ssGSEA score", add = "none", 
          palette = c("High" = mycol[1], "Low"=mycol[2])
) +
  coord_flip() +    
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5) + 
  theme(legend.position = "top", 
        legend.key.size = unit(10,'pt'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black'),
        axis.text.y = element_text(color = "black", hjust = 1, size = 10, lineheight = 2),
        plot.title = element_text(hjust = 0, size = 10),
        text = element_text(family = "")) +
  theme(axis.title.x = element_text(size =10, lineheight = 2),
        axis.title.y = element_text(size =10, lineheight = 2)) +
  annotate(geom = 'text', label = stat_multi$P.sig,hjust = 0, size = unit(2.5,'mm'), color = 'black',
           x = unique(data2$Pathway), y=max(data2$Score)-0.3)  
ggsave("TIDEother.pdf", height = 4, width = 4.5)
