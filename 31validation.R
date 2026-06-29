read.countData <- read.table("bindGEO_GPL96_remove.txt",header=T,sep="\t",comment.char="",check.names=F)
rt2 =as.matrix(read.countData)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
countData2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)

group=read.csv("group.csv",header=T,sep=",",comment.char="",check.names=T,row.names = 1)
conNum=length(group[group=='Normal'])   
treatNum=length(group[group=='PTC'])     

gene="NOS3"   
data=t(countData2[gene,,drop=F])
data=data[rownames(group),]
exp=cbind(data, group)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
NT=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=NT)
comp=combn(NT,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


library(ggplot2)
library(ggpubr)
library(gghalves)
ggplot(exp, aes(x = Type, y = gene, color = Type, fill = Type)) +
  geom_boxplot(outlier.shape = 2, width = 0.3, notch = F, size = 0.7) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = alpha(mycol,0)) +
  labs(y = NULL, x = gene) +
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  stat_compare_means(comparisons = list(c("Normal","PTC")),
                     method = "wilcox.test", paired = F,  
                     label = "p.format")

ggsave("NOS3.pdf",height = 3.5,width = 2)





