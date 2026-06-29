library(GSVA)
library(limma)
library(GSEABase)


immuneScore=function(expFile=null, gmtFile=null, project=null){
  rt=read.table(expFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  mat=avereps(mat)
  mat=mat[rowMeans(mat)>0,]
  
  geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
  
  param=ssgseaParam(mat,geneSet)
  ssgseaScore=gsva(param, verbose = TRUE)
  normalize=function(x){
    return((x-min(x))/(max(x)-min(x)))}
 
  ssgseaOut=normalize(ssgseaScore)
  ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
  write.table(ssgseaOut, file=paste0(project, ".score.txt"), sep="\t", quote=F, col.names=F)
}

immuneScore(expFile="total.normalize.txt", 
            gmtFile="19.allTCpathway", 
            project="ssgsea")


ssgsea = read.table("ssgsea.score.txt", header=T, sep="\t", check.names=F, row.names=1)
ssgsea = t(ssgsea)

library(reshape2)
cluster=read.table("1.3clusterall_log.txt", header=T, sep="\t", check.names=F, row.names=1)  #Cluster数据04，或者Group数据12
rt2 = cbind(ssgsea, cluster)
data2 = melt(rt2, id.vars = "Cluster")
colnames(data2) = c("Cluster","Pathway","Score")

data2$Pathway = sort(data2$Pathway, decreasing = F)  


kruskal_data <- data.frame(Pathway = colnames(rt2)[1:ncol(rt2)-1])
for (i in 1:(ncol(rt2)-1))  
  {print(i)
  kruskal_data[i,2] <- kruskal.test(rt2[,i] ~ Cluster, data = rt2)[["p.value"]]}
kruskal_data$fdr <- p.adjust(kruskal_data$V2, method = "fdr")
kruskal_data$BH <- p.adjust(kruskal_data$V2, method = "BH")
colnames(kruskal_data) <- c("Pathway","p.value","adj.p.value")   
kruskal_data$p.value3 = ifelse(kruskal_data$p.value<0.001,"<0.001",sprintf("%.03f", kruskal_data$p.value))
kruskal_data$p.signif <- ifelse(kruskal_data$adj.p.value < 0.001, "***",
                                ifelse(kruskal_data$adj.p.value < 0.01,"**",
                                       ifelse(kruskal_data$adj.p.value < 0.05, "*", "NS")))
kruskal_data$p.sci = ifelse(kruskal_data$p.value>0.001,sprintf("%.03f", kruskal_data$p.value),
                            format(kruskal_data$p.value, scientific = TRUE, digits = 1))


library(tidyr)
data_long <- rt2 %>% 
  pivot_longer(cols = colnames(rt2)[1:27], names_to = "Pathway", values_to = "Score")
anova_results <- data.frame(Pathway = character(), P_value = numeric(),stringsAsFactors = FALSE)
for (p in unique(data_long$Pathway)) {
  pathway_data <- data_long %>% filter(data_long$Pathway == p)
  anova_result <- aov(data_long$Score ~ data_long$Cluster, data = pathway_data)
  p_value <- summary(anova_result)[[1]]$"Pr(>F)"[1]
  anova_results <- rbind(anova_results, data.frame(Pathway = p, P_value = p_value))
}
print(anova_results)

data2$Cluster <- factor(data2$Cluster, levels = c("C1","C2",'C3')) 


library(ggridges)
library(ggplot2)
ggplot(data2, aes(y = Pathway, x = Score, fill = Cluster, color = Cluster)) +
  geom_density_ridges(scale = 1, alpha = 0.7, size = 0.1,
                      jittered_points = F , rel_min_height = 0.01,
                    
                      position = position_points_jitter(height = 0.2, width = 0.1), 
                  
  ) +
  scale_fill_manual(values = c("C1"=alpha(mycol[1],0.8),"C2"=alpha(mycol[2],0.5),"C3"=alpha(mycol[3],0.5))) +
  scale_color_manual(values = c("C1"= alpha(mycol[1],0.8),"C2"=alpha(mycol[2],0.5),"C3"=alpha(mycol[3],0.5)))+
  scale_y_discrete(limits = rev(levels(data2$Pathway))) +
  xlim(0,1.1) + 
  labs(title = "", x = "Score", y = "") +
  theme_ridges() +
  theme(legend.position = "top") +
  annotate(geom = 'text',label = paste0(kruskal_data$p.sci, kruskal_data$p.signif),, hjust = 0, size = unit(3,'mm'), color = 'black', #fontface = faceCell,
           y = unique(data2$Pathway), x=1.01)  


ggsave("Cluster.ssGSEA.pdf",height = 9, width = 10)

