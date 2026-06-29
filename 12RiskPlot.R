library(ggrisk)
library(rms)    
library(pheatmap)
library(survival)
library(ggplot2)

inputFile = "trainRisk.txt"
rt=read.table(inputFile, header=T, sep="\t", check.names=F,row.names=1)    
rt=rt[order(rt$RiskScore),]      
rt$patient_id = seq(nrow(rt))  
riskClass=rt[,"Risk"]
lowLength=length(riskClass[riskClass=="Low"])
highLength=length(riskClass[riskClass=="High"])
lowMax=max(rt$RiskScore[riskClass=="Low"])
rt$RiskScore2 <- rt$RiskScore
rt$RiskScore2[rt$RiskScore2 > 5] = 5
rt$color <- c(rep("#589cd6", lowLength), rep("#E57259", nrow(rt) - lowLength))
p1 <- ggplot(rt, aes(x = rt$patient_id, y = rt$RiskScore2)) +
  geom_point(aes(color = I(color)), shape = 20, size = 1.5) +
  scale_color_manual(values = c("#589cd6", "#E57259"), labels = c("Low-risk", "High-risk")) +
  labs(x = "Patients (increasing Score)", y = "Score") +
  theme_bw() +
 geom_hline(yintercept = lowMax, linetype = "dashed") +
  geom_vline(xintercept = lowLength, linetype = "dashed") + 
  theme(legend.position = "right", legend.box.spacing = unit(0.2, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10))+
  guides(color = guide_legend(title = "Group", override.aes = list(shape = 20, size=2.5)
  )) +  theme(axis.title = element_text(size=10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 10,face="bold")) +
  ggtitle("Entire cohort")    

p2 <- ggplot(rt, aes(x = rt$patient_id)) +
  geom_point(aes(y = rt$DSS, color = as.factor(rt$Event)),shape = 20, size = 2.5,alpha = 1) +
  scale_color_manual(values = c("#589cd6","#E57259"), labels = c("Censored","Event")) +
  labs(x = "Patient (increasing Score)", y = "Follow-up time (year)") +
  geom_vline(xintercept = lowLength, linetype = "dashed") +
  theme_bw() +
  guides(color = guide_legend(title = "Status",override.aes = list(shape = 20, size=2.5))) +
  theme(legend.position = "right", legend.box.spacing = unit(0.2, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10),
        axis.title = element_text(size=10)) 

library(patchwork)
p1 + p2 + 
  plot_layout(nrow = 2,   
              heights = c(1, 1),  
              guides = "keep") &
  theme(legend.position = 'right',
        legend.box.spacing = unit(0.1,"cm"),
        legend.spacing = unit(0,'cm'),
        legend.justification = "centre",
        legend.key.size = unit(0.4,"cm")) 

ggsave("train.pdf",height = 4,width = 6)  

