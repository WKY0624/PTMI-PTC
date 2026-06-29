library(TCGAbiolinks)

project <- getGDCprojects()$project_id 
project <- project[grep("TCGA-", project)]
project <- "TGCA-THCA"

url = "http://biocc.hrbmu.edu.cn/TIP/RCodeAndData/pancancerData/"

dir.create ("TIP")

dir.create ("TIP/Immune activity scores/") 
dir.create("TIP/Immune cell infiltration/") 
dir.create("TIP/SignatureGenes.Expression/")

for(proj in project){
  message(proj)
  cancer = unlist(strsplit(proj, "-"))[2] 
  download.file(url = paste0(url,cancer,"/ssGSEA.normalized.score.txt"),
                destfile = paste0("TIP/Immune activity scores/", proj, ".txt"))
  download.file(url = paste0(url,cancer,"/CIBER_",cancer,"_lm14_allsample_Result.txt"),
                destfile = paste0("TIP/Immune cell infiltration/", proj, ".txt"))
  download.file(url = paste0(url,cancer,"/SignatureGenes.Expression.txt"),
                destfile = paste0("TIP/SignatureGenes.Expression/",proj,".txt"))
}

library(ggpubr)
rt=read.table("TGCA-THCA.txt",sep="\t",header=T,row.names=1,check.names=F)   
group=sapply(strsplit(colnames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
colnames(rt)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(rt))

ClusterFile="1.3clusterall_log.txt" 
Cluster = read.table(ClusterFile,sep="\t",check.names=F,row.names=1,header=T)
Cluster = Cluster[,ncol(Cluster),drop=F]
RiskFile="totalRisk.txt"
Risk = read.table(RiskFile,sep="\t",check.names=F,row.names=1,header=T)
Risk = Risk[,c(ncol(Risk)-1,ncol(Risk)),drop=F]

samSample = intersect(colnames(rt),row.names(Cluster))
data = rt[,samSample,drop=F]
Cluster = Cluster[samSample,,drop=F]
merge=cbind(t(data),Cluster,Risk) 


merge$Step4 <- rowMeans(merge[,4:20])
merge <- merge[,-c(4:20)]
merge <- merge[, c("Step1", "Step2", "Step3", "Step4", "Step5", "Step6", "Step7", "Cluster","RiskScore")]
merge <- merge[order(merge$RiskScore, decreasing = F),]

C1 <- merge[merge$Cluster=="C1",,drop=F]
C1 <- C1[order(C1$RiskScore, decreasing = T),]
C2 <- merge[merge$Cluster=="C2",,drop=F]
C2 <- C2[order(C2$RiskScore, decreasing = T),]
C3 <- merge[merge$Cluster=="C3",,drop=F]
C3 <- C3[order(C3$RiskScore, decreasing = T),]

merge.order <- rbind(C1,C2,C3)

data_plot <- merge.order[,-c(8,9)]
data_plot <- as.data.frame(t(data_plot))

outTab = data.frame()
pSig = c()
pNum = c()
for(i in colnames(merge.order[,1:(ncol(merge.order)-2)])){
  rt1 = merge.order[,c(i,"Cluster")]
  colnames(rt1)= c("expression","Cluster")
  ksTest = kruskal.test(expression ~ Cluster, data = rt1)
  pValue = ksTest$p.value
  p1 <- ifelse(pValue<0.001,"***",
               ifelse(pValue<0.01,"**",
                      ifelse(pValue<0.05,"*","ns"))) 
  p2 <- ifelse(pValue>0.001,round(pValue,3),format(pValue, scientific = TRUE, digits = 3))
  if(pValue < 1){        
    outTab=rbind(outTab,cbind(rt1,gene=i))
    print(pValue)}
  pSig[i] <- p1
  pNum[i] <- p2
}

pSig = as.data.frame(pSig)
pNum = as.data.frame(pNum)
pSig$stepSig = paste0(pNum$pNum,pSig$pSig)

rownames(data_plot) <- pSig$stepSig

annotation_col = data.frame(merge.order$Cluster)
colnames(annotation_col) = "Cluster"
rownames(annotation_col) <- rownames(merge.order)

annotation_row = data.frame(factor(c(1:7)))
colnames(annotation_row) = "Cycle"
rownames(annotation_row) = rownames(data_plot)
head(annotation_row)


library(RColorBrewer)

color_scheme <- "Purples"
col7 <- brewer.pal(7, color_scheme)

ann_colors = list(
  Cluster = c(C1 = mycol[1], C2 = mycol[2],C3=mycol[3]),
  Cycle = c("1" = col7[1], "2" = col7[2], "3" = col7[3], "4" = col7[4], "5" = col7[5], "6" = col7[6], "7" = col7[7]))

gap1 = nrow(C1)
gap2 = nrow(C2)
gap3 = nrow(C3)



library(ComplexHeatmap)
pdf("TIPheatmap.pdf", width = 6.5, height = 3)
pheatmap(data_plot,   
         cluster_cols = F,
         cluster_rows = F,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#4474c4",1.5), "white", rep("#EB7E60",1.5)))(100),
         fontsize_col = 8,
         fontsize_row = 10,
         show_colnames = F,
         gaps_col = c(gap1, gap1+gap2),
         cellwidth = 0.5, 
         cellheight = 20)

dev.off()



step4 = merge[,grep("Step4|Cluster|Risk", colnames(merge))]
outTab = data.frame()
pSig = c()
pNum = c()
for(i in colnames(step4[,1:(ncol(step4)-3)])){
  rt1 = step4[,c(i,"Cluster")]
  colnames(rt1)= c("expression","Cluster")
  Test = kruskal.test(expression ~ Cluster, data = rt1) 
  pValue = Test$p.value
  p1 <- ifelse(pValue<0.001,"***",
               ifelse(pValue<0.01,"**",
                      ifelse(pValue<0.05,"*","ns"))) 
  p2 <- ifelse(pValue>0.001,round(pValue,3),format(pValue, scientific = TRUE, digits = 1))
  if(pValue < 1){        #⭐️
    outTab=rbind(outTab,cbind(rt1,step=i))
    print(pValue)}
  pSig[i] <- p1
  pNum[i] <- p2
}

pSig = as.data.frame(pSig)
pSig$stepSig = paste0(rownames(pSig),pSig$pSig)
pNum = as.data.frame(pNum)

outTab = data.frame()
pSig = c()
pNum = c()
for(i in colnames(step4[,1:(ncol(step4)-3)])){
  rt1 = step4[,c(i,"Risk")]
  colnames(rt1)= c("expression","Risk")
  Test = wilcox.test(expression ~ Risk, data = rt1)
  pValue = Test$p.value
  p1 <- ifelse(pValue<0.001,"***",
               ifelse(pValue<0.01,"**",
                      ifelse(pValue<0.05,"*","ns"))) 
  p2 <- ifelse(pValue>0.001,round(pValue,3),format(pValue, scientific = TRUE, digits = 1))
  if(pValue < 1){        #⭐️
    outTab=rbind(outTab,cbind(rt1,step=i))
    print(pValue)}
  pSig[i] <- p1
  pNum[i] <- p2
}

pSig = as.data.frame(pSig)
pSig$stepSig = paste0(rownames(pSig),pSig$pSig)
pNum = as.data.frame(pNum)


library(ggradar)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)

data <- outTab
table(data$step)

new_data <- data %>%
  group_by(step, Cluster) %>%
  mutate(Means = mean(expression),
         Median = median(expression))

plot_data <- unique(new_data[,c(2,3,5)]) %>%     
  pivot_wider(names_from = step, values_from = Median) %>%
  column_to_rownames("Cluster")

data2 <- plot_data
data2$group <- rownames(data2)

data3 <- data2 %>% 
  as_tibble(rownames = "Cluster") %>% 
  select(1:18)    

colnames(data3) = c("Cluster",pSig$stepSig)
colnames(data3) <- gsub(".recruiting", "", colnames(data3))  
colnames(data3) <- gsub("Step4.", "", colnames(data3))  

ggradar(data3,
        font.radar = 'sans',
        base.size = 2,
        
           grid.min = min(plot_data),
        grid.mid = (max(plot_data) + min(plot_data))/2,
        grid.max = max(plot_data),
        
      
        label.gridline.min = F, 
        label.gridline.mid = F,
        label.gridline.max = F,
        grid.label.size	= 3,    
        grid.line.width	= 0.7,  
        gridline.min.linetype = 2,   
        gridline.mid.linetype	= 2, 
        gridline.max.linetype = 2, 
        gridline.min.colour	= alpha('grey',0.7),
        gridline.mid.colour	= alpha('grey',0.7),
        gridline.max.colour	= alpha('grey',0.7),
        
        axis.line.colour = alpha('grey',0.7), 
        axis.label.offset	= 1.1,   
        axis.label.size	= 3,   
        
        background.circle.colour = 'gray',
        background.circle.transparency	= 0.07,
        
        group.colours = c("#EB7E60","#aad09d","#8FB4DC"),   
        group.point.size = 2,      
        group.line.width = 0.8,     
        fill = F,   
        fill.alpha = 0.2,   
        
        legend.title = "Cluster",
        legend.text.size = 1,
        legend.position = "top") +
  
  theme_void() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.direction = 'horizontal',
        legend.position = 'top',
  )


ggsave('TIP.Step4.Cluster.Radar.pdf', width = 5, height = 4)

