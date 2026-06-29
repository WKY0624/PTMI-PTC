
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

dir="Training Data"
dir(dir)


myexp = read.table("total.normalize.txt", sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = F)  
myexp = as.matrix(myexp)
myexp = log2(myexp+1)
dim(myexp)

group = read.table("totalRisk.txt", sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F)  
group = group[,ncol(group),drop=F]
head(group)

exp = readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
exp[1:4,1:4]
dim(exp)   
drug = readRDS(file = file.path(dir,"CTRP2_Res.rds"))
drug <- exp(drug) 
ggboxplot(melt(drug[ , 1:4]), "Var2", "value") 
drug[1:4,1:4]
dim(drug)  
identical(rownames(drug),colnames(exp))

calcPhenotype(trainingExprData = exp,
              trainingPtype = drug,
              testExprData = myexp,
              batchCorrect = 'eb',    
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'homogenizeData',
 )


library(data.table)
testPtype <- read.csv("DrugPredictions.csv", row.names = 1,check.names = F)
testPtype[1:4, 1:4]
dim(testPtype)
identical(colnames(testPtype),colnames(drug)) 
library(stringr)
a = t(rbind(drug,testPtype))
a = a[,c(1,546,2,547,3,548,4,549)]
par(mfrow = c(2,2))
plot(a[,1],a[,2])
plot(a[,3],a[,4])
plot(a[,5],a[,6])
plot(a[,7],a[,8])


ggboxplot(reshape2::melt(drug[,1:4]), "Var1", "value") 
round(apply(drug[1:4 ,], 1, function(x){
  return(c(
    head(sort(x)),
    tail(sort(x))
  ))
}),2)


apply(drug[ 1:4 ,], 1, function(x){ 
  names(x)=gsub('_[0-9]*','',colnames(drug))
  return(c(
    names(head(sort(x))),
    names(tail(sort(x)))
  ))
})


rm(list = ls())
library(tidyr)
library(dplyr)
library(ggplot2)

Drugfile = "DrugPredictions.csv"
Riskfile = "totalRisk.txt"  #12
RAIfile = "27.RAI.txt"

group = read.table(Riskfile, sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
group = group[,ncol(group),drop=F]

RAI = read.table(RAIfile, sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F )
sameSample = intersect(row.names(group), row.names(RAI))
RAIgroup = RAI[sameSample,,drop=F]
group = group[sameSample,,drop=F]
RAIgroup = RAIgroup[RAIgroup$RAI %in% c("Sensitive","Refractory"),,drop=F]  
RAIgroup$SampleID = rownames(RAIgroup)

drug <- read.csv(Drugfile, row.names = 1, check.names = F)   
sameSample2 = intersect(row.names(RAIgroup), row.names(drug))
drug = drug[sameSample2,,drop=F]
group = group[sameSample2,,drop=F]


sameSample=intersect(row.names(group), row.names(drug))
group = group[sameSample,"Risk",drop=F]
group$SampleID = rownames(group)
drug = drug[sameSample,,drop=F]
combined.data <- cbind(data.frame(SampleID = rownames(drug)), as.data.frame(drug))
combined.data <- pivot_longer(combined.data, -SampleID, names_to = "Drug", values_to = "IC50")
combined.data <- left_join(combined.data, group, by = "SampleID")


drug.data <- combined.data %>% filter(grepl("axitinib", Drug, ignore.case = TRUE))

re_outlier<-function(x,na.rm = TRUE,...){
  qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...)  
  h<-1.5*IQR(x,na.rm = na.rm)
  y<-x
  y[x<(qnt[1]-h)]<-NA
  y[x>(qnt[2]+h)]<-NA
  y}
drug.data  <- drug.data %>%
  group_by(Risk)%>%
  mutate(IC50 = re_outlier(IC50))
drug.data <- drug.data[complete.cases(drug.data),]

drug.diff <- drug.data %>%
  group_by(Risk) %>%
  summarize(Mean = mean(IC50), SD = sd(IC50),
            Median = median(IC50), Q1 = quantile(IC50, 0.25), Q3 = quantile(IC50, 0.75), Count = n())
print(drug.diff)

stat.data <- wilcox.test(drug.data$IC50[drug.data$Risk == "Low"],
                         drug.data$IC50[drug.data$Risk == "High"],
                         paired = F, alternative = 'two.sided')

library(ggplot2)

ggplot(drug.data, aes(x = Risk, y = IC50, fill = Risk))+
  scale_fill_manual(values = mycol) + 
  geom_violin(aes(color = Risk), width = 0.5, size = 1, alpha=0.1) +
  scale_color_manual(values = mycol)+  
  geom_boxplot(notch = F, outlier.size = -1, width = 0.5,lwd=0.5, color="black", alpha = 0.6)+ 
  ylab("IC50") +
  xlab("Group")  +
  theme_classic() +
  theme(
    axis.ticks.y = element_line(size=0.5, color="#333333"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.y = unit(0.3,"cm"),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_blank()) +
  labs(tag = paste0(unique(drug.data$Drug), "\nP ", ifelse(stat.data$p.value<0.001, paste0("= ", format(stat.data$p.value, scientific = TRUE, digits = 3)), paste0("= ",round(stat.data$p.value,3)))))+
  theme(plot.tag.position = c(0.55,0.9),
        plot.tag = element_text(size = 10, color = 'black',
                                vjust = 0,   
                                hjust = 0.5,  
                                lineheight = 1.48   
        )) 

unique(drug.data$Drug)
ggsave('axitinib.pdf', height = 3.2, width = 3.5)


results <- data.frame(Drug = character(0), p_value = numeric(0),
                      Median.High = numeric(0), Median.Low = numeric(0),
                      Mean.High = numeric(0), Mean.Low = numeric(0))

for (drug.name in unique(combined.data$Drug)) {
  drug.data <- combined.data %>% filter(Drug == drug.name)
    re_outlier<-function(x,na.rm = TRUE,...){
    qnt<-quantile(x,probs = c(0.25,0.75),na.rm = na.rm,...) 
    h<-1.5*IQR(x,na.rm = na.rm)
    y<-x
    y[x<(qnt[1]-h)]<-NA
    y[x>(qnt[2]+h)]<-NA
    y}
  drug.data  <- drug.data %>%
    group_by(Risk)%>%
    mutate(IC50 = re_outlier(IC50))
  drug.data <- drug.data[complete.cases(drug.data),]
    drug.diff <- drug.data %>%
    group_by(Risk) %>%
    summarize(Mean = mean(IC50), SD = sd(IC50),
              Median = median(IC50), Q1 = quantile(IC50, 0.25), Q3 = quantile(IC50, 0.75), Count = n())
  
  drug.low <- drug.data$IC50[drug.data$Risk == "Low"]
  drug.high <- drug.data$IC50[drug.data$Risk == "High"]
    wilcox.result <- wilcox.test(drug.high,drug.low)
  
  
  results <- rbind(results, data.frame(Drug = drug.name, 
                                       Median.High = drug.diff$Median[drug.diff$Risk == "High"],
                                       Median.Low = drug.diff$Median[drug.diff$Risk == "Low"],
                                       Mean.High = drug.diff$Mean[drug.diff$Risk == "High"],
                                       Mean.Low = drug.diff$Mean[drug.diff$Risk == "Low"],
                                       p_value = wilcox.result$p.value))
}

results <- results[order(results$p_value), ]
results <- results[order(results$Median.High), ]

top_10_drugs <- head(results, 10)
top_30_drugs <- head(results, 30)



score = read.table(Riskfile, sep = '\t', header = T, row.names = 1, check.names = F, stringsAsFactors = F)
score = score[,ncol(score)-1,drop=F]
score = score[sameSample,,drop=F]
score$ID = rownames(score)
drug2 <- drug
drug2$ID = rownames(drug2)
drug.score = merge(score, drug2, by.x = "ID", by.y = "ID")
rownames(drug.score) = drug.score$ID
drug.score = subset(drug.score,select = -c(ID)) 
score1 = score[,ncol(score)-1,drop=F]

Spearman=data.frame()
for(agents in colnames(drug)){
  for(RiskScore in colnames(score1)){
    x=as.numeric(drug[,agents])
    y=as.numeric(score1[,RiskScore])
    corT=cor.test(x,y,method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    Spearman = rbind(Spearman,cbind(Agents = agents, cor, text, pvalue))
  }
}

data = merge(Spearman, results, by.x = "Agents", by.y = "Drug")
data <- data[order(data$p_value), ]
data30 <- head(data, 30)   
data10 <- head(data, 10)  

data30$Agents = factor(data30$Agents, levels = data30$Agents)
data30$cor = as.numeric(data30$cor)
data30$p_value = as.numeric(data30$p_value)
data30$size <- with(data30, 2 + (5 - 2) * ((Mean.Low - min(Mean.Low)) / (max(Mean.Low) - min(Mean.Low))))

data30 = data30[order(data30$p_value,decreasing = T), ]
data30$log10P = as.numeric(-log10(data30$p_value))


data10$Agents = factor(data10$Agents, levels = data10$Agents)
data10$cor = as.numeric(data10$cor)
data10$p_value = as.numeric(data10$p_value)
data10$size <- with(data10, 2 + (5 - 2) * ((Mean.Low - min(Mean.Low)) / (max(Mean.Low) - min(Mean.Low))))
data10 = data10[order(data10$p_value,decreasing = T), ]
data10$log10P = as.numeric(-log10(data10$p_value))

data30_long <- data30 %>%
  pivot_longer(c(Mean.High,Mean.Low),
               names_to = "group",
               values_to = "IC50")
data10_long <- data10 %>%
  pivot_longer(c(Mean.High,Mean.Low),
               names_to = "group",
               values_to = "IC50")

data30_long_mm <- data30_long %>%
  group_by(Agents) %>%
  mutate(x_min = min(IC50),
         x_max = max(IC50))
data10_long_mm <- data10_long %>%
  group_by(Agents) %>%
  mutate(x_min = min(IC50),
         x_max = max(IC50))

ggplot(data30, aes(x=reorder(Agents, log10P), y=cor)) +
  geom_segment(aes(x=reorder(Agents, log10P), xend=Agents, y=0, yend=cor,  
                   #linetype=factor(data$tumortype,levels = c("Primary","Metastatic"))
  ), 
  alpha = 1, linewidth = 0.3
  ) +
  geom_point(aes(size = Mean.Low, color = log10P),pch = 20) +
  geom_hline(yintercept = 0, linetype = 2, color = 'gray20', linewidth = 0.3) +
  scale_size(range = c(2,5)) +
  scale_colour_gradient2(low = alpha('#d2d4f5',1), mid = alpha('#fdbfca',1), high = alpha("#d87070",1),   #"dodgerblue2"  "#fb3e35"
                        
                         limits = c(min(data30$log10P), max(data30$log10P)),
                         breaks = c(5,6,7,8,9,10),
                         midpoint = 5,
                         name = "-log10(pvalue)") +  
  ylim(-0.5,0.5) +
  coord_flip() +
  theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(axis.text = element_text(size = 8, color = '#333333'),
        axis.text.x = element_text(vjust = 0, hjust = 0.5),  
  ) +
  xlab("") +
  ylab("Spearman correlation") +
  theme(legend.position = 'right',
        legend.key.size = unit(10,'pt'),
        legend.title = element_text(size = 8, lineheight = 4),
        legend.text = element_text(size = 8, lineheight = 5),)  +
  guides(size = guide_legend(title = "IC50", order = 0, keyheight=1)) 


ggsave('top30.pdf', width = 5, height = 5)

