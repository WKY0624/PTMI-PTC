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
  
  param=ssgseaParam(mat,
                    geneSet,
                    assay = NA_character_,
                    annotation = NULL,
                    minSize = 1,
                    maxSize = Inf,
                    alpha = 0.25,
                    normalize = TRUE
                    )
  ssgseaScore=gsva(param)
  
  normalize=function(x){
    return((x-min(x))/(max(x)-min(x)))}
  ssgseaOut=normalize(ssgseaScore)
  ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
  write.table(ssgseaOut, file=paste0(project, ".score.txt"), sep="\t", quote=F, col.names=F)
}

immuneScore(expFile="total.normalize.txt", 
            gmtFile="immune.gmt", 
            project="total")


library(limma)
library(reshape2)
library(ggpubr)


data=read.table("total.score.txt", header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

risk=read.table("totalRisk.txt", header=T, sep="\t", check.names=F, row.names=1)  #12

sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
rt=cbind(data,risk[,c("RiskScore","Risk")]) 
rt=rt[,-(ncol(rt)-1)]


immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
          "Neutrophils","pDCs","T_helper_cells",
          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
rt1=rt[,c("Risk",immCell)]   
data=melt(rt1,id.vars=c("Risk"))   
colnames(data)=c("Risk","Type","Score")   
data$Risk=factor(data$Risk, levels=c("High","Low"))

data$Type = sort(data$Type, decreasing = T)   
library(rstatix)  
library(ggpubr)
wilcox_test(aDCs~Risk, data = rt1, alternative = 'two.sided')
wilcox_data <- data.frame(Genesymbol=colnames(rt1)[2:ncol(rt1)])
for (i in 2:ncol(rt1)){
  print(i)
  wilcox_data[i-1,2] <- wilcox.test(rt1[,i] ~ Risk, data = rt1, 
                                  alternative = 'two.sided', 
                                  exact = FALSE)[["p.value"]]
}

wilcox_data$fdr <- p.adjust(wilcox_data$V2, method = "fdr")

colnames(wilcox_data) <- c("Genesymbol","p.value","adj.p.value")   
wilcox_data$p.value3 = ifelse(wilcox_data$p.value<0.001,"<0.001",sprintf("%.03f", wilcox_data$p.value))
wilcox_data$p.signif = ifelse(wilcox_data$adj.p.value<0.001,"***",
                              ifelse(wilcox_data$adj.p.value<0.01,"**",
                                     ifelse(wilcox_data$adj.p.value<0.05,"*","")))

Cell = unique(data$Type)
immCell2=paste0(immCell,"(",wilcox_data$p.value,")")  
colorsCell = ifelse(wilcox_data$p.value<0.05,"#fb3e35","black")
faceCell = ifelse(wilcox_data$p.value<0.05,2,1)

p1<-ggboxplot(data, x='Type', y= "Score", color = "Risk",
          notch = F, size = 0.4, width = 0.6, outlier.shape = 3,
          xlab="Immune cell",ylab="ssGSEA score", add = "none", 
          palette = c("High" = mycol[1], "Low"=mycol[2])
) +
  coord_flip() +   
  scale_y_continuous(limits = c(0,1.)) +
  rotate_x_text(50) +
 
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5) + 
  theme(legend.position = "top", 
        legend.key.size = unit(10,'pt'),
      
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
       
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black'),
        axis.text.y = element_text(face = rev(faceCell),
                                  
                                   hjust = 1, 
                                   size = 10,
                                   lineheight = 2),
        plot.title = element_text(hjust = 0, size = 10),
        text = element_text(family = "")) +
  theme(axis.title.x = element_text(size =10, lineheight = 2),
        axis.title.y = element_text(size =10, lineheight = 2)) +
  annotate(geom = 'text', label = wilcox_data$p.signif,hjust = 0, size = unit(3,'mm'), color = 'black',fontface = faceCell,
           x = unique(data$Type), y=1)  
stat_compare_means(aes(group=`Risk`),
                   symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                    symbols = c("***", "**", "*", "")),
                   method="wilcox.test",  
                   label.y = 1.05,
                   label.x = ,
                   label.y.npc = 'centre',
                   label = "p.format") 

ggsave('Cellrisk.pdf',height = 4, width = 5)



immFunc=c("APC_co_inhibition","APC_co_stimulation","CCR",
          "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
          "MHC_class_I","Parainflammation","T_cell_co-inhibition",
          "T_cell_co-stimulation","Type_I_IFN_Reponse")
rt2=rt[,c("Risk",immFunc)]
data2=melt(rt2,id.vars=c("Risk"))
colnames(data2)=c("Risk","Type","Score")
data2$risk=factor(data2$Risk, levels=c("High","Low"))

data2$Type = sort(data2$Type, decreasing = T)   
library(rstatix)  
library(ggpubr)
wilcox_test(APC_co_inhibition ~ Risk, data = rt2)
wilcox_data2 <- data.frame(Genesymbol=colnames(rt2)[2:ncol(rt2)])
for (i in 2:ncol(rt2)){
  print(i)
  wilcox_data2[i-1,2] <- wilcox.test(rt2[,i] ~ Risk, data = rt2, exact = FALSE)[["p.value"]]
}

wilcox_data2$fdr <- p.adjust(wilcox_data2$V2,method = "fdr")
data_control2 <- rt2[c(which(rt2$Risk == "Low")),]
data_trt2 <- rt2[c(which(rt2$Risk == "High")),]
wilcox_data2$foldchange <- colMeans(data_trt2[2:ncol(data_trt2)])-colMeans(data_control2[2:ncol(data_control2)])

colnames(wilcox_data2) <- c("Genesymbol","p.value","adj.p.value","logFC")
wilcox_data2$p.value3 = ifelse(wilcox_data2$p.value<0.001,"<0.001",sprintf("%.03f", wilcox_data2$p.value))
wilcox_data2$p.signif = ifelse(wilcox_data2$adj.p.value<0.001,"***",
                               ifelse(wilcox_data2$adj.p.value<0.01,"**",
                                      ifelse(wilcox_data2$adj.p.value<0.05,"*","")))


immFunc2=paste0(immFunc,"(",wilcox_data2$p.value,")")  
colorsFunc = ifelse(wilcox_data2$p.value<0.05,"#fb3e35","black")
faceFunc = ifelse(wilcox_data2$p.value<0.05,2,1)
p2 <- ggboxplot(data2, x='Type', y= "Score", color = "Risk",
                notch = F, size = 0.4, width = 0.6, outlier.shape = 3,
                xlab="Immune function",ylab="ssGSEA score", add = "none",
                palette = c("High" = mycol[1], "Low"=mycol[2])
) +
  coord_flip() +   
  scale_x_discrete(labels = sort(immFunc,decreasing = T),     
                   position = ) +   
  scale_y_continuous(limits = c(0,1)) +
  rotate_x_text(50) +
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5) + 
  theme(legend.position = "top", 
        legend.key.size = unit(10,'pt'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black'),
        axis.text.y = element_text(
          face = rev(faceFunc),
          hjust = 1, 
          size = 10,
          lineheight = 2),
        plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(family = "")) +
  theme(axis.title.x = element_text(size =10, lineheight = 2),
        axis.title.y = element_text(size =10, lineheight = 2)) +
  annotate(geom = 'text', label = wilcox_data2$p.signif,hjust = 0, size = unit(3,'mm'), color = 'black',fontface = faceFunc,
           x = unique(data2$Type), y=1)   
stat_compare_means(aes(group=`Risk`),
                   symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                    symbols = c("***", "**", "*", "")),
                   method="wilcox.test",              
                   label.y = 1.05,
                   label.x = ,
                   label.y.npc = 'centre',
                   label = "p.format") 


ggsave('Functionrisk.pdf',height = 4, width = 5.1)

