library(limma)
library(reshape2)
library(ggplot2)
library(dplyr)
library(grid)

expFile="total.normalize.txt"   
riskFile="totalRisk.txt"
geneFile = "26.ImmunomodulatorsList1.txt"   
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), rownames(data))
data1=t(data[sameGene,])

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data1), row.names(risk))
data1=data1[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F]

outTab=data.frame()
for(checkpiont in colnames(data1)){   
  for(gene in colnames(risk)){        
    x=as.numeric(data1[,checkpiont])    
    y=as.numeric(risk[,gene])         
    corT=cor.test(x,y,method="spearman")   
    cor=corT$estimate
    pvalue=corT$p.value
    text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
    outTab=rbind(outTab,cbind(Gene=gene, checkpiont=checkpiont, cor, text, pvalue))
  }
}



outTab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$checkpiont=factor(outTab$checkpiont, levels=colnames(data1))
outTab$cor=as.numeric(outTab$cor)
outTab$pvalue=as.numeric(outTab$pvalue)
ysort = unique(outTab$checkpiont)
xface = c('plain','plain','plain','plain','plain','plain','plain','plain','plain','plain','plain','plain','bold')


library(aplot)
library(tidyr)
p <- ggplot(outTab, aes(Gene, checkpiont)) + 
  geom_tile(aes(fill = cor), size = 0)+   
  scale_fill_gradient2(low = mycol2[2] , mid = "white", high = mycol2[1],  
                       breaks=c(-0.5,0,0.5), labels = c(-0.5,0,0.5), 
                       limits=c(-0.56,0.72),
                       name = "rho") + 
  geom_text(aes(label=text),col ="#333333",size = 3) +   
  theme_void() +    
  scale_y_discrete(limits=factor(sort(ysort, decreasing = T)),position = 'left')+
  theme_test(base_size = 10, base_line_size = 0.3, base_rect_size = 0.5) +
  theme(axis.title = element_blank(), 
        axis.ticks.x = element_line(linewidth = 0.3),  
        axis.ticks.y = element_line(linewidth = 0.3), 
        axis.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = xface, color = "#333333", margin = margin(0.2,0,0,0, 'cm')),   #x轴字体
  ) +      
  theme(legend.position = "bottom", legend.key.height = unit(2,"mm"), legend.key.width = unit(10,"mm"),
        legend.direction = "horizontal") +
  scale_x_discrete(position = "bottom") +        
  guides(fill='none')   

table(geneRT$Immunomodulators)
geneRT$Gene <- factor(geneRT$Gene, levels = sort(unique(geneRT$Gene), decreasing = F))
left <- geneRT$Gene %>%  as.data.frame() %>% 
  mutate(group = rep(unique(geneRT$Immunomodulators), times = c(46,21,24,24,33,22,46))) %>%    
  mutate(p="")  %>%   
  ggplot(aes(x=p, y=geneRT$Gene, fill = group)) + 
  geom_tile() +
  scale_fill_manual(values = mycol[1:7]) +
  scale_y_discrete(position="right") +  
  theme_void() +
  theme(axis.text.y = element_text(colour = '#333333', size= 10, hjust = 1, margin = margin(0,0.1,0,0, 'cm'))) +
  theme(legend.position = "none") +
  xlab(NULL) + ylab(NULL) +
  labs(fill='group')

p %>% insert_left(left, width = .08) 


ggsave('ImmunomodulatorsCor.pdf', width =3.6, height = 13)  #mod1（137）:3.5x20，mod2（78）：3.5x13\11\12


dt = data.frame(Gene = c('A','B','C','D','E','F','G'), Geneset = c('Chemokines','HLA','Inhibitor','Interferons','Interleukins','Other cytokines','Stimulator'), cor = 1:7)  
for_legend <- ggplot(dt, aes(Gene, Geneset)) + 
  geom_tile(aes(fill = cor), colour = "white", size = 0.5)+
  scale_fill_gradient2(low = mycol2[2], mid = "white", high = mycol[1],  
                       breaks=c(-0.5,0,0.5),labels=NULL,   
                       limits=c(-1,1),
                       name = "rho") +
  guides(fill = guide_colorbar(
    title.position = 'left',
    title.theme = element_text(size = 8,face = "bold",colour = "black"),
    label = T,  
    label.theme = element_text(size = 8,face = "plain",colour = "black"),
    raster = T,
    frame.colour = NULL,    
    barwidth = unit(20,"mm"),
    barheight = unit(3,"mm"),
    nbin = 50,  
    ticks = T,   
    draw.ulim = T,  
    draw.llim = T,  
  )) +
  theme(legend.direction = "horizontal")

legend1 <- cowplot::get_legend(for_legend)
grid.newpage()
grid.draw(legend1)


library(ComplexHeatmap)
legend2 <- Legend(labels = c('Chemokines','HLA','Inhibitor','Interferons','Interleukins','Other cytokines','Stimulator'), 
                  title = "Immunomodulators", title_gap = unit(3,'mm'),
                  row_gap = unit(2,"mm"),
                  legend_gp = gpar(fill = color))
draw(legend2)

