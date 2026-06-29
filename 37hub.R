read.countData <- read.table("04.Counts_vst.txt",header=T,sep="\t",comment.char="",check.names=F)
rt2 =as.matrix(read.countData)
rownames(rt2)=rt2[,1]
exp2=rt2[,2:ncol(rt2)]
dimnames2=list(rownames(exp2),colnames(exp2))
countData2=matrix(as.numeric(as.matrix(exp2)),nrow=nrow(exp2),dimnames=dimnames2)

gene="TYMS"   
data=t(countData2[gene,,drop=F])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))


cli=read.table("06.clinical477.txt",header=T,sep="\t",comment.char="",check.names=F,row.names = 1)
cli=cli[rownames(data),]
data=cbind(data,cli)
data=as.data.frame(data)
exp=data[,grep('TYMS|Age|Gender|CombinedDisease|Tcategory2|Tcategory|Ncategory|Mcategory|HistologicalType|ETE|Multifocality|TumorDiameter|BRAF|RAS|TERT|RAIResponse|NewTumorEvent',colnames(data))]
exp$Tcategory<- factor(exp$Tcategory, levels = c("T1",'T2','T3','T4'))
exp <- exp[!is.na(exp$Tcategory), ] 
group = levels(exp$Tcategory)
my_comparisons = list(group)


library(ggplot2)
library(ggpubr)
p <- ggplot(exp,aes(x = Tcategory, y = TYMS, color = Tcategory))+   
  geom_violin(   
    alpha = 0.8,    
    scale = 'width',  
    trim = TRUE)+
  geom_boxplot(mapping=aes(x=Tcategory,y=TYMS,colour=Tcategory,fill=Tcategory),               
               alpha = 0.5,                 
               size=1.5,                 
               width = 0.3)+   
  geom_jitter(mapping=aes(x=Tcategory,y=TYMS,colour = Tcategory),             
              alpha = 0.3,size=3)+  
  scale_fill_manual(limits=c("T1",'T2','T3','T4'), 
                    values =mycol)+  
  scale_color_manual(limits=c("T1", 'T2','T3','T4'),                    
                     values=mycol)+   
  geom_signif(mapping=aes(x=Tcategory,y=TYMS),           
              comparisons = list(                              
              c("T1", "T2"),                                 
              c("T1", "T3"),
              c('T1','T4'),
              c('T2','T3'),
              c('T3','T4')),   
            
              map_signif_level=F,              
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0),               
              y_position = c(9.4,9.6,9.5,9.4),              
              size=1,            
              textsize = 4,              
              test = "t.test",
              color='#5b679b')+ 
  theme_bw()+
  labs(x="Tcategory",y="TYMS") 
p   
ggsave("Tcategory.pdf",height = 5,width = 9)


ggplot(exp, aes(x = Age, y = TYMS, color = Age, fill = Age)) +
  geom_boxplot(outlier.shape = 8, width = 0.3, notch = T, size = 0.7) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = alpha(mycol,0.3)) +
  labs(y = NULL, x ='Age')+
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.ticks.x = element_line(color = "black")
  )+
  stat_compare_means(comparisons = list(c("<55",">=55")),
                     method = "t.test", paired = F,  
                     label = "p.format")


p <- ggplot(exp, aes(x = Age, y = TYMS, fill = Age)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9,
              alpha = 0.7, color = NA) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4, shape = 3, color = 'white') +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label = "p.signif", bracket.size = 0.6, tip.length = 0.03,
                     size = 5) +
  scale_fill_manual(values = mycol) +
  labs(x = "Age group", y = "TYMS expression") +
  theme_classic(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "None",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

