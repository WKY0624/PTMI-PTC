library(ggpubr) 
library(ggplot2)
library(gghalves)
library(reshape) 
library(rstatix)  
tciaFile="29.TCIA-THCA_230310_tidy.tsv"  
riskFile="totalRisk.txt" 

ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)
Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(ips), row.names(Risk))
ips=ips[sameSample, , drop=F]
Risk=Risk[sameSample, "Risk", drop=F]
data=cbind(Risk,ips)
high  = data[data$Risk=='High',]
low = data[data$Risk=='Low',] 

ks.test(scale(data$ips_ctla4_neg_pd1_pos),'pnorm')  
ks.test(scale(low$ips_ctla4_pos_pd1_neg),'pnorm')  


wilcox.test(ips_ctla4_neg_pd1_pos ~ Risk, data = data, alternative = 'two.side')
stat <- data.frame(TCIA=colnames(data)[2:ncol(data)])
for (i in 2:ncol(data)){
  print(i)
  stat[i-1,2] <- wilcox.test(data[,i] ~ Risk, data = data, 
                             alternative = 'two.side',
                             exact = FALSE)[["p.value"]]}

stat$fdr <- p.adjust(stat$V2, method = "fdr")
stat$bon <- p.adjust(stat$V2, method = 'bonferroni')

colnames(stat) <- c("TCIA","P.value","FDR","Bonferroni")
data_long = melt(data, id.vars = ("Risk"))
colnames(data_long) = c("RiskGroup", "TCIA", "Score")

data_long <- data_long %>% group_by(TCIA)
data_long$Score <- as.numeric(data_long$Score)

stat.data <- pairwise_wilcox_test(data = data_long, Score~RiskGroup, paired = F, alternative = 'two.side') %>% 
  add_xy_position(x = 'TCIA')
stat.data$p.scient <- format(stat.data$p.adj, scientific = TRUE)
stat.data$p.round3 <- round(stat.data$p.adj,3)
stat.data$Pvalue <- ifelse(stat.data$p<0.001, sprintf(stat.data$p.scient), sprintf("%.03f", stat.data$p))

ggplot()+
  geom_half_violin(data = data_long %>% filter(RiskGroup == "High"),
                   aes(x = TCIA, y = Score), side = "l",size= 0,
                   colour="white", fill=mycol[1], alpha = 0.2, width = 1,nudge = 0.01,  
                   position = position_dodge(width = 0)
  ) +
  geom_half_violin(data = data_long %>% filter(RiskGroup == "Low"),
                   aes(x = TCIA,y = Score), side = "r", size= 0,
                   colour="white", fill=mycol[2], alpha = 0.2, width = 1.1, nudge = 0.01, 
                   position = position_dodge(width = 0)
  ) +
  geom_half_boxplot(data = data_long %>% filter(RiskGroup == "High"),
                    aes(x = TCIA, y = Score), width = 0.2, 
                    colour=mycol[1], #lwd= 0.4, 
                    outlier.shape = NA, outlier.size = 0.8, outlier.stroke = T,
                    fill=mycol[1],side = "l", alpha = 0.6, nudge = 0.01, errorbar.draw = F,
                    position = position_dodge(width = 1))+
  geom_half_boxplot(data = data_long %>% filter(RiskGroup == "Low"),
                    aes(x = TCIA, y = Score), width = 0.2, 
                    colour=mycol[2], #lwd= 0.4, 
                    outlier.shape = NA, outlier.size = 0.8, outlier.stroke = T,
                    fill=mycol[2],side = "r", alpha = 0.6, nudge = 0.01, errorbar.draw = F,
                    position = position_dodge(width = 1)
  ) +
 
  geom_line(data = data_long, aes(x = TCIA, y = Score, color = RiskGroup),
            stat = 'summary', fun=median, lty =1, size =2,
            position = position_dodge(width = 0.1)) +
  scale_color_manual(values = mycol,name="Group",labels = c("High","Low")) +
  stat_pvalue_manual(
    stat.data, 
    y.position = c(10.5,10.5,10.5,10.5),
    label = 'P = {Pvalue}\n{p.adj.signif}',  
    bracket.size = 0.3,
    bracket.shorten = 0.15,  
    tip.length = 0.01, 
  ) +
  xlab("") +
  ylab("TCIA Score") +
  ylim(4,12) + 
  scale_x_discrete(labels = c('IPS','PD1 blocker','CTLA4 blocker','PD1+CTLA4 blocker'))+
  theme_classic(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5)+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black'),    
    axis.ticks =element_line(linewidth = 0.3),
    axis.text.y = element_text(color = 'black', hjust = 1, # 左对齐
                               size = 10, lineheight = 1),
    plot.title = element_text(hjust = 0.5, size = 10),
    text = element_text(family = ""),
    legend.position = "top",  
    legend.key.size = unit(10,'pt'), 
    legend.justification = "centre")

ggsave('TCIA.pdf', height = 4.5, width = 6)
