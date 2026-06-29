score=read.table("TMEscores_TPM.txt", header=T, sep="\t", check.names=F, row.names=1)
score=as.matrix(score)
row.names(score)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(score))
score=avereps(score)
score=score[,1:4]  

read.countData <- read.table("04.Counts_vst.txt",header=T,sep="\t",comment.char="",check.names=F)
rt2 =as.matrix(read.countData)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
countData2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)

gene="TYMS"   
data=t(countData2[gene,,drop=F])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))

data=data[rownames(score),'TYMS']


data = cbind(data,score)
colnames(data)[colnames(data) == "data"] <- "TYMS"
data=as.data.frame(data)

SpearmanR <- cor(data$TYMS, data$StromalScore,method="spearman",use="complete.obs")
SpearmanP <- cor.test(data$TYMS, data$StromalScore, method="spearman", use="complete.obs")

ggplot(data, aes(x=TYMS, y=StromalScore)) +
  geom_point(color = mycol2[1], alpha=0.7, pch=20, size=2) +
  geom_smooth(method=lm , formula = y ~ x, 
             
              color=mycol2[1], fill=mycol2[1], alpha = 0.3, se=TRUE) +
  theme_test(base_line_size = 0.3)+
  ylab("StromalScore") + 
  xlab("TYMS") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 10),
    axis.title.y = element_text(vjust = 0),
    axis.text = element_text(color = "gray30",size = 9)
  ) +   
  xlim(min(data$TYMS), max(data$TYMS)) + 
  ylim(min(data$StromalScore), max(data$StromalScore))+  
  annotate("text", x = max(data$TYMS), y = max(data$StromalScore)-500, fontface = 1, hjust = 1, label = paste0("rho = ",round(SpearmanR,3), "\n","P = ", format(SpearmanP$p.value,scientific = TRUE)))
ggsave("StromalScore.pdf", height = 3, width = 3.2)

