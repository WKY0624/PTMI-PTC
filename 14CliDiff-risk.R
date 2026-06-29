risk=read.table("totalRisk.txt",header=T, sep="\t", check.names=F, row.names = 1)
clin = read.table("06.clinical477.txt",header=T, sep="\t", check.names=F, row.names = 1)
merge = merge(risk, clin, by = "row.names", all = TRUE)
rownames(merge) = merge[,1]
merge = merge[,-1]

rt1 = merge[,c("RiskScore","Age","Gender","RadiationHistory","CombinedDisease")]

rt2 = merge[,c("RiskScore","Tcategory2","Ncategory","Mcategory","HistologicalType","ETE","Multifocality")]
rt2$ETE = ifelse(rt2$ETE == "None","None","ETE")


rt3 = merge[,c("RiskScore","BRAF","RAS","TERT","GeneFusion")]
rt3$BRAF = ifelse(rt3$BRAF == "Mutation","BRAF","non-BRAF")
rt3$RAS = ifelse(rt3$RAS == "Mutation","RAS","non-RAS")
rt3$TERT = ifelse(rt3$TERT == "Mutation","TERT","non-TERT")


rt4 = merge[,c("RiskScore","RAIResponse","NewTumorEvent","TNM","ATA")]
rt4$NewTumorEvent = ifelse(rt4$NewTumorEvent == "None","None","DSS")

library(tidyverse)
long_data <- pivot_longer(rt4, cols = c(RAIResponse,NewTumorEvent,TNM,ATA), names_to = "class", values_to = "class_temp")

result_data <- cbind(long_data$RiskScore, long_data$class_temp)  

colnames(result_data) <- c("Riskscore", "class")
result_data = as.data.frame(result_data)
result_data$Riskscore = as.numeric(result_data$Riskscore)


library(ggplot2)
library(ggpubr)   
custom_order = c("Not received","Sensitive","Refractory","None","DSS","Stage I","Stage II","Stage III~IV","Low","Intermediate","High")



result_data$class = factor(result_data$class, levels = custom_order)
result_data= na.omit(result_data)    


library(ggbeeswarm)
ggplot(result_data, aes(x = class, y=Riskscore, color = class, fill = class)) + 
  geom_beeswarm(cex = 1, size = 2 , pch = 16, alpha = 0.4, stroke = 0, na.rm = TRUE,
                priority = "none", 
  ) + 
  geom_boxplot(outlier.shape = NA, width = 0.3,  color="#333333", size = 0.4) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = alpha(mycol,0)) +
  
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.4) +
  
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        axis.text.x = element_text(size = 10, colour = 'black',angle = 45, hjust = 1),
        legend.position = "none",
  ) +
  
  labs(y = "Riskscore", x = "") +  
  ylim(1,1.4) +
  
  stat_compare_means(
    method = "wilcox.test", 
    method.args = list(alternative = "two.side"),
    comparisons = list(
      c("Not received","Sensitive"),c("Not received","Refractory"),c("Sensitive","Refractory"),
      c("None","DSS"),c("Stage I","Stage II"),c("Stage I","Stage III~IV"),c("Stage II","Stage III~IV"),
      c("Low","Intermediate"),c("Low","High"),c("Intermediate","High")
    ),
    label = "p.signif", 
    hide.ns = F,
    step.increase = F)

ggsave("RiskScore4.pdf",height = 4, width = 7)
