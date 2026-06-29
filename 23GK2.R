library("clusterProfiler")
library("org.Hs.eg.db")  
library("enrichplot")
library("ggplot2")
library("stringr")

pvalueFilter=0.05     
qvalueFilter=0.05           

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
rt=read.table("M1_diff_1.3MRGsall_stat.xls", header=T, sep="\t", check.names=F)   

gene=read.table("AI_1.3Genes.txt", header=T, sep="\t", check.names=F)
gene=gene$x
rt <- rt[rt$ID %in% gene, ]

genes=as.vector(rt[,1])
entrezIDs=mget(genes, 
               org.Hs.egSYMBOL2EG, 
               ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
genes=entrezIDs[entrezIDs!="NA"]        

kk1=enrichGO(gene = genes,
            OrgDb = org.Hs.eg.db,   
            pvalueCutoff =0.05,  
            qvalueCutoff = 0.05,   
            ont="all",   
            readable =T)
GO=as.data.frame(kk1)

write.table(GO,file="GO.total.txt",sep="\t",quote=F,row.names = F)


kk2 <- enrichKEGG(gene = genes, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05,minGSSize = 5)
KEGG=as.data.frame(kk2)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)


showNum=20
if(nrow(GO)<showNum){
  showNum=nrow(GO)
}
showNum=20
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}


pdf(file="bubble.pdf",width = 10,height =7)
bub=dotplot(kk1,
            x = "GeneRatio",
            color = colorSel,   
            size = "Count", 
            showCategory = showNum,
            orderBy = "GeneRatio", 
            label_format = 100,  
            split="ONTOLOGY",
            font.size = 10) + facet_grid(ONTOLOGY~., scale='free')  
print(bub)
dev.off()

pdf(file="KEGG.bub.pdf",width = 9,height = 15)
dotplot(kk2, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()


pdf(file="barplot.pdf",width = 10,height =7)
bar=barplot(kk1, 
            drop = TRUE, 
            showCategory =showNum,
            split="ONTOLOGY",
            color = colorSel,
            label_format = 100, 
            font.size = 10) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()


pdf(file="KEGG.bar.pdf",width = 9,height = 13)
barplot(kk2, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()
library('enrichplot')
pdf(file="KEGG.for.pdf",width = 13,height = 10)
cnetplot( kk2,
          label_format = 25,
          showCategory = 25,
          node_label = "category",    
          colorEdge = F,
          color.params = list(
            edge = "#F9B797",
            gene = "#9271B1", 
            category = c("low" = "#FAE2C1", "high" = "#E4908B") 
          ))
dev.off()


library(tidyverse)

go_top<-as.data.frame(GO)%>%
  group_by(ONTOLOGY)%>%
  slice_head(n=5)%>%
  arrange(desc(qvalue))%>%
  ungroup()%>%
  dplyr::select(ONTOLOGY,everything())

kegg_top<-as.data.frame(kk2)%>%
  dplyr::arrange(desc(qvalue))%>%
  dplyr::slice(1:7)%>%
  dplyr::select(ID:Count)%>%
  dplyr::mutate(ONTOLOGY="KEGG")%>%
  dplyr::select(ONTOLOGY,everything())

data<-rbind(go_top,kegg_top)%>%
  dplyr::mutate(ONTOLOGY=factor(ONTOLOGY,levels=rev(c("BP","CC","MF","KEGG")),ordered=TRUE))%>%
  dplyr::arrange(desc(ONTOLOGY),-log10(qvalue))


data$geneID<-gsub("/",".",data$geneID)
data$CountNumber<-data$Count/1e3


plot<-ggplot(data)+
  geom_bar(aes(x=-log10(qvalue),y=interaction(Description,ONTOLOGY),
               fill=ONTOLOGY),stat="identity")+
  scale_fill_manual(values=color,name="ONTOLOGY")+
  geom_text(aes(x=0.1,y=interaction(Description,ONTOLOGY),
                label=Description),size=3,hjust=0,color="black")+
  geom_text(aes(x=0.1,y=interaction(Description,ONTOLOGY),
                label=geneID),size=2,hjust=0,vjust=2.5,color="black")+
  geom_point(aes(x=-max(Count)/20,y=interaction(Description,ONTOLOGY),
                 size=Count,fill=ONTOLOGY),shape=21)+
  geom_text(aes(x=-max(Count)/20,y=interaction(Description,ONTOLOGY),label=Count),size=3)+
  scale_size(range=c(4,8),guide=guide_legend(override.aes=list(fill="black")))+
  guides(fill=guide_legend(reverse=TRUE))+
  labs(x="-log10(FDR)",y="Description")+
  theme(
    legend.title=element_text(color="#000000",size=12),
    legend.text=element_text(color="#000000",size=10),
    axis.text.x=element_text(color="#000000",size=12),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_text(color="#000000",size=14),
    panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    legend.background=element_blank()
  )
ggsave(plot=plot,filename="GOKEGG.pdf",height=7,width=7)


