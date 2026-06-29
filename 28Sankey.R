library(ggalluvial)
library(ggplot2)
library(dplyr)


clusterFile = "1.3clusterall_log.txt"    
scoreFile="totalRisk.txt" 
clinicalFile = "06.clinical477.txt"

clusterRT = read.table(clusterFile, header = T, sep = '\t', check.names = F)
scoreRT = read.table(scoreFile, header = T, sep = '\t', check.names = F)
clinicalRT = read.table(clinicalFile, header = T, sep = '\t', check.names = F)
TCfi = read.table("06.time_DSS.txt", header = T, sep = '\t', check.names = F)
clinicalRT=merge(clinicalRT,TCfi,by.x = 'Sample', by.y = "ID")


scoreRT$new_column<-rownames(scoreRT)
colnames(scoreRT)[which(colnames(scoreRT) == "new_column")] <- "id"

score = scoreRT[,c('id','Risk'),drop=F]
clinical = clinicalRT[,c('Sample','Event','RAIResponse'),drop = F]
TIDE = read.table("TIDE.csv", header = T, sep = ',', check.names = F)
TIDE=TIDE[,c('Patient','Responder'),drop = F]

merge2 = merge(score, clinical, by.x = 'id', by.y = "Sample")  
merge3 = merge(clusterRT,merge2, by.x = 'ID', by.y = "id")
merge4 = merge(merge3,TIDE,by.x='ID',by.y='Patient')

data = merge4[,c(2:6),drop=F]
data=data[,-3]
data = na.omit(data) 

data$Cluster = ifelse(data$Cluster == 'C1','C1',
                      ifelse(data$Cluster == 'C2','C2','C3'))
data$Risk=ifelse(data$Risk=='High','High-risk','Low-risk')
data$Responder=ifelse(data$Responder=='True','Respond','Non-respond')

colnames(data)[colnames(data)=="Cluster"] <- "Cluster"   
colnames(data)[colnames(data)=="Risk"] <- "Group"   
colnames(data)[colnames(data)=="RAIResponse"] <- "RAI"  
colnames(data)[colnames(data)=="Responder"] <- "TIDE"  


corLodes = to_lodes_form(data, axes = 1:ncol(data), id = "Cohort")

corLodes$stratum = factor(corLodes$stratum, 
                          levels = c('C1','C2','C3','High-risk','Low-risk',
                                     'Refractory','Sensitive','Not received','Non-respond','Respond'))

ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 1)) + 
  geom_flow(width = 1/10, aes.flow = "forward") +   
  geom_stratum(alpha = 1, width = 3.5/10,color = 'white',linetype=1, lwd=0.3) +  
  scale_fill_manual(values = mycol) +
  geom_text(stat = "stratum", size = 3.3, color="black") + 
  xlab("") + ylab("") + 
  theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank(),
        axis.text.x = element_text(color = 'black', size = 10)) +
  theme(panel.grid =element_blank()) +   
  theme(panel.border = element_blank()) +   
  theme(legend.position = "none") + 
  ggtitle("") + 
  guides(fill = FALSE)  


ggsave('Sankey.pdf', height = 3.5, width = 4)
