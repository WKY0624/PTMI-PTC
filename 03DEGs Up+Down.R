library(ggplot2)
library(scales)
library(dplyr)

data=read.table("M1_diff_1.3MRGsall_stat.xls",header=T,sep="\t",comment.char="",check.names=F)

my_colors <- colorRampPalette(c("#FFC5C5", "#FF8C8C", "#FF4C4C"))
scale_fill_gradientn(colours = my_colors(3), name = "adj.P.Val")

ggplot(data, 
       aes(x=reorder(ID, log2FoldChange), 
           y=log2FoldChange, fill=padj)) +
  geom_bar(stat="identity") +  
  coord_flip() +  
  scale_fill_gradient(low="#7d3293", high="#D8BFD8", name="adj.P.Val") +  
  theme_classic() +  
  labs(    title = "Gene Expression Log Fold Change",    
           subtitle = "Adjusted P-Values Indicated by Fill Color",    
           x = "Gene",    
           y = "logFC"  ) + 
  theme(axis.text.y = element_text(size=8),  
        axis.line = element_line(colour = "black", size = 0.5),  
        axis.line.x = element_line(colour = "black", size = 0.5),
      legend.position = "top",  
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
      plot.subtitle = element_text(hjust = 0.5, size = 12))+
  geom_text(aes(label = sprintf("%.2f", log2FoldChange)),
            hjust = 1.1, vjust = 0, size = 3, colour = "black") +  
  theme(panel.grid.major.x = element_line(colour = "grey90", linetype = "dashed")) 


data$log10PValue <- -log10(data$pvalue)

data$group <- case_when(data$log2FoldChange > 1 & data$pvalue < 0.05 ~ "Up",
                      data$log2FoldChange < -1 & data$pvalue < 0.05 ~ "Down",
                      abs(data$log2FoldChange) <= 1 ~ "None",
                      data$pvalue >= 0.05 ~ "None")
up<- filter(data,group=="Up") %>% top_n(10,log2FoldChange)
down<- filter(data,group=="Down") %>% top_n(10,abs(log2FoldChange))


up <- up %>% arrange(desc(log2FoldChange))
up$ID <- factor(up$ID,levels = rev(up$ID),ordered = T)
down <- down %>% arrange(log2FoldChange)
down$ID <- factor(down$ID,levels = rev(down$ID),ordered = T)
df<- rbind.data.frame(up,down)
df$ID <- factor(df$ID,levels = rev(df$ID),ordered = T)


p1<- ggplot(df, aes(x = log2FoldChange,y = ID,fill = log10PValue)) +
  geom_col(color = "white",width=0.85,linewidth=0.7) +
  labs(x="log2FC",y="",fill="-log10 Pvalue")
p1

lx <- c(rep(-0.2,10),rep(0.2,10))

top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
mytheme <- theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        legend.title.position = "left",
        legend.key.height=unit(0.8, "cm"),,
        legend.key.width=unit(0.4, "cm"),,
        legend.title = element_text(size = 9,hjust = 0.5, angle = 90),
        axis.line = element_line(linewidth = 0.5,lineend = "square"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),units="inches"))

p2 <- p1+ geom_text(aes(label = ID,x = lx,color = group,
                        hjust = ifelse(log2FoldChange < 0, 0, 1)),
                    show.legend = F,size = 3,fontface = "bold")+
  scale_color_manual(values = mycolor)+
  scale_fill_gradient(low = "#ffffff",high = "#ff8c00",
                      limits = c(0,100),breaks = c(0,20,40,60,80,100))
p2+mytheme
p3<- p2+ scale_x_continuous(breaks = seq(-4, 8, by = 2),
                             limits= c(-4, 8),expand = c(0,0)) 
p3+mytheme
p4 <- p3 + geom_text(aes(label = round(log2FoldChange,2),
                         hjust = ifelse(log2FoldChange < 0, -0.2, 1.2),vjust = 0.5),
                     color="#ffffff",size = 3,fontface = "bold")+
  geom_vline(xintercept = 0, color = "#ffffff", lwd = 0.5)
p4+mytheme


mytheme1<- theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        legend.title.position = "left",
        legend.key.height=unit(0.5, "cm"),,
        legend.key.width=unit(0.4, "cm"),,
        legend.title = element_text(size = 9,hjust = 0.5, angle = 90),
        plot.margin=unit(x=c(0.5,0.2,0,0.2),units="inches"))

f1<- ggplot(df, aes(x = log2FoldChange,y = ID,fill = log10PValue)) +
  geom_col(data=up,color = "white",width=0.9,linewidth=0.7) +
  labs(x="log2FC",y="",fill="-log10 Pvalue")+
  geom_text(data=up,aes(label = ID,x = lx[1:10],color = group,
                        hjust= ifelse(log2FoldChange < 0, 0, 1)),
            show.legend = F,size = 3,fontface = "bold",)+
  scale_color_manual(values = mycolor[2])+
  scale_fill_gradient(low = "#ffffff",high = mycolor[2],
                      limits= c(0,120),breaks = c(0,30,60,90,120))+
  scale_x_continuous(breaks = seq(-4, 8, by = 2),
                     limits= c(-4, 8),expand = c(0,0)) +
  geom_text(data=up,aes(label = round(log2FoldChange,2),
                        hjust= ifelse(log2FoldChange < 0, -0.2, 1.2),vjust = 0.5),
            color="#ffffff",size = 3,fontface = "bold")+
  mytheme1
f1
mytheme2<- theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_blank(),
        legend.title.position = "left",
        legend.key.height=unit(0.5, "cm"),,
        legend.key.width=unit(0.4, "cm"),,
        legend.title = element_text(size = 9,hjust = 0.5, angle = 90),
        axis.line = element_line(linewidth = 0.5,lineend = "square"),
        plot.margin=unit(x=c(0,0.2,0.1,0.2),units="inches"))

f2<- ggplot(df, aes(x = log2FoldChange,y = ID,fill = log10PValue)) +
  geom_col(data=down,color = "white",width=0.9,linewidth=0.7) +
  labs(x="log2FC",y="",fill="-log10 Pvalue")+
  geom_text(data=down,aes(label = ID,x = lx[11:20],color = group,
                          hjust= ifelse(log2FoldChange < 0, 0, 1)),show.legend = F,
            size= 3,fontface = "bold")+
  scale_color_manual(values = mycolor[1])+
  scale_fill_gradient(low = "#ffffff",high = mycolor[1],
                      limits= c(0,40),breaks = c(0,10,20,30,40))+
  scale_x_continuous(breaks = seq(-4, 8, by = 2),
                     limits= c(-4, 8),expand = c(0,0)) +
  geom_text(data=down,aes(label = round(log2FoldChange,2),
                          hjust= ifelse(log2FoldChange < 0, -0.2, 1.2),vjust = 0.5),
            color="#ffffff",size = 3,fontface = "bold")+
  mytheme2
f2
library(cowplot)
plot_grid(f1, f2, ncol = 1, align = "v")




