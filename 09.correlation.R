
library(corrr)
library(dplyr)

rt=read.table("AI_1.3GenesExp.xls",header=T,sep="\t",comment.char="",check.names=F,row.names = 1)
rt=rt[-c(13:15),]
rt=t(rt)

rt %>%
  correlate() %>%
  rearrange() %>%
  network_plot(colours = c('skyblue1', "white","indianred2"))


library(corrplot)
library(ggplot2)
library(ggcorrplot)
library(vcd)
library(ggrepel)
library(psych)

data=rt

dim(data)
data<-as.matrix(data) 
data=data.frame(scale(data))
head(data)

cor_matrix <- corr.test(data,method = "spearman")

round(data, 2)
r <- cor_matrix$r 
p <- cor_matrix$p 
pdf('Correlation.pdf',height=7,width=7.5)

corrplot(r, method="circle",          
         title = "Spearman",    
         type="lower",      
         col=colorRampPalette(c("#4474c4",'#8FB4DC','#e99c93','indianred2'))(50),,          
         outline = T,        
         diag = TRUE,  
         mar = c(0,0,1,0),        
         bg="white",   
         add = FALSE,       
         is.corr = TRUE,        
         addgrid.col = "darkgray",         
         addCoef.col = NULL,        
         addCoefasPercent = FALSE,      
         order = "original",     
         hclust.method = "complete",   
         addrect = NULL,        
         rect.col = "black", 
         rect.lwd = 2, 
         tl.pos = "l", 
         sig.level = 0.05,            
         insig = "label_sig",           
         pch.cex = 1.8)
corrplot(p, title = "",                 
         method = "number",                
         outline = F,      
         add = TRUE,          
         type = "upper",        
         order="original",         
         col=colorRampPalette(c("#4474c4",'#8FB4DC','gray80','#e99c93','indianred2'))(50),         
         diag=FALSE,    
         tl.pos="n",       
         cl.pos='n' 
        )
dev.off()






