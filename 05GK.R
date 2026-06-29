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

genes=as.vector(rt[,1])
up_genes <- as.vector(rt$ID[rt$log2FoldChange > 0 & rt$padj< 0.05]) 
down_genes <- as.vector(rt$ID[rt$log2FoldChange < 0 & rt$padj< 0.05])
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
write.table(GO,file="GO.totalup.txt",sep="\t",quote=F,row.names = F)
write.table(GO,file="GO.totaldown.txt",sep="\t",quote=F,row.names = F)


kk2 <- enrichKEGG(gene = genes, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05,minGSSize = 5)
KEGG=as.data.frame(kk2)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)
write.table(KEGG,file="KEGGup.txt",sep="\t",quote=F,row.names = F)
write.table(KEGG,file="KEGGdown.txt",sep="\t",quote=F,row.names = F)


showNum=20
if(nrow(GO)<showNum){
  showNum=nrow(GO)
}
showNum=10
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

GOtotal<-read.table("GO.total.txt", header=T, sep="\t", check.names=F)
GOup<-read.table("GO.totalup.txt", header=T, sep="\t", check.names=F)
GOdown<-read.table("GO.totaldown.txt", header=T, sep="\t", check.names=F)

GOtotal$new_column<-"TotalSet"
colnames(GOtotal)[which(colnames(GOtotal) == "new_column")] <- "Group"
rownames(GOtotal)=GOtotal[,2]


GOup$new_column<-"UP"
colnames(GOup)[which(colnames(GOup) == "new_column")] <- "Group"

GOdown$new_column<-"DOWN"
colnames(GOdown)[which(colnames(GOdown) == "new_column")] <- "Group"

GO_GK<-rbind(GOup,GOdown)
GO_GK <- GO_GK%>%
  dplyr::select(1,3,4,11,12,13,14) 

KEGGtotal<-read.table("KEGG.txt", header=T, sep="\t", check.names=F)
KEGGup<-read.table("KEGGup.txt", header=T, sep="\t", check.names=F)
KEGGdown<-read.table("KEGGdown.txt", header=T, sep="\t", check.names=F)

KEGGtotal$new_column<-"TotalSet"
colnames(KEGGtotal)[which(colnames(KEGGtotal) == "new_column")] <- "Group"

KEGGup$new_column<-"UP"
colnames(KEGGup)[which(colnames(KEGGup) == "new_column")] <- "Group"

KEGGdown$new_column<-"DOWN"
colnames(KEGGdown)[which(colnames(KEGGdown) == "new_column")] <- "Group"
KEGG_GK<-rbind(KEGGup,KEGGdown)
KEGG_GK<- KEGG_GK%>%
  dplyr::select(4,5,12,13,14,15) 
KEGG_GK$new_column<-"KEGG"
colnames(KEGG_GK)[which(colnames(KEGG_GK) == "new_column")] <- "ONTOLOGY"


GK<-rbind(GO_GK,KEGG_GK)
GK <- GK %>%
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::slice_max(GeneRatio, n = 10) %>%
  dplyr::ungroup()

GK<-GK%>%
  dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("BP","CC","MF","KEGG"), ordered = T))%>%
  dplyr::mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))

p=ggplot(GK,aes(x=Group,
                y=Description,
                size=Count,
                colour=qvalue))+
  geom_point(shape=16)+
  labs(x="Group",y="Pathway")+
  scale_colour_continuous(
    name="Enrichment",
    low = "#E9CDDF", high = "#ce1256")+
  scale_radius(
    range=c(5,9),
    name="Size")+
  guides(
    color=guide_colorbar(order=1),
    size=guide_legend(order=2)
  )+theme_bw()
p
ggsave(filename = "udbarplot.pdf",
       height = 12,
       width = 8)


library(tidyverse)
library(ggfun)
library(patchwork)
Up_result<-read.table("GO.totalup.txt", header=T, sep="\t", check.names=F)%>%
  dplyr::select(1,2,3,13)
Down_result<-read.table("GO.totaldown.txt", header=T, sep="\t", check.names=F)%>%
  dplyr::select(1,2,3,13)


plot_data <- inner_join(Up_result,
                        Down_result,
                        by = c("ONTOLOGY", "ID", "Description")) %>%
  dplyr::rename(Up = Count.x,
                Down = Count.y) %>%
  dplyr::group_by(ONTOLOGY) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::ungroup()

GO <- c("BP","MF","CC")
GO2 <- c("Biological Process","Molecular Function","Cellular Component")

plot_list_out <- lapply(1:3, function(x){
  p <- plot_data %>%
    dplyr::filter(ONTOLOGY == GO[x]) %>%
    dplyr::mutate(sum = Up + Down) %>%
    dplyr::arrange(sum) %>%
    dplyr::mutate(Description = factor(Description, levels = Description, ordered = T)) %>%
    dplyr::select(-sum) %>%
    tidyr::pivot_longer(cols = c(Up, Down), names_to = "Kind", values_to = "Number") %>%
    ggplot(aes(x = Number, y = Description)) + 
    geom_bar(aes(fill = Kind), stat = "identity", width = 0.75) + 
    labs(x = "", y = GO2[x]) + 
    scale_fill_manual(values = c("Up" = "#F3B082",
                                 "Down" = "#9DC3E6")) + 
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) + 
    theme_bw() +
    theme(
      axis.text = element_text(size = 10, color = "#000000"),
      axis.title = element_text(size = 15, color = "#000000"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x  = element_line(linetype = 2, linewidth = 1),
      panel.grid.minor.x  = element_blank(),
      panel.border = element_rect(linewidth = 1),
      legend.background = element_roundrect(color = "#969696")
    )
  
  return(p)
  
})

p_combine <- plot_list_out[[1]]/plot_list_out[[2]]/plot_list_out[[3]] + 
  plot_layout(guides = "collect",
              heights = c(2,0.8,0.5))


ggsave(filename = "GO.pdf",
       plot = p_combine,
       height = 10,
       width = 8)


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

#合并
data<-rbind(go_top,kegg_top)%>%
  dplyr::mutate(ONTOLOGY=factor(ONTOLOGY,levels=rev(c("BP","CC","MF","KEGG")),ordered=TRUE))%>%
  dplyr::arrange(desc(ONTOLOGY),-log10(qvalue))

data<-go_top%>%
  dplyr::mutate(ONTOLOGY=factor(ONTOLOGY,levels=rev(c("BP","CC","MF")),ordered=TRUE))%>%
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


data <- as.data.frame(GO) %>%
  dplyr::select('Description', 'ONTOLOGY', 'GeneRatio', 'qvalue') %>%
  mutate(len = str_length(Description)) %>%
  dplyr::filter(len < 120) %>%
  rowwise() %>%
  mutate(GeneRatio = round(eval(parse(text = GeneRatio)), 3) * 100) %>%
  arrange(qvalue) %>%
  group_by(ONTOLOGY) %>%
  mutate(ID = 1:n()) %>%
  top_n(10, wt = -ID) %>% {
    tmp <- as.data.frame(KEGG) %>%
      dplyr::select('Description', 'GeneRatio', 'qvalue') %>%
      rowwise() %>%
      mutate(GeneRatio = round(eval(parse(text = GeneRatio)), 3) * 100,
             ONTOLOGY = 'KEGG') %>%
      arrange(qvalue) %>%
      group_by(ONTOLOGY) %>%
      mutate(ID = 1:n()) %>%
      top_n(10, wt = -ID)
    rbind(., tmp)
  } %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF", "KEGG")),
         Description = str_wrap(Description, width = 60),
         Description = factor(Description, levels = rev(Description)))

ggplot(data, aes(-log10(qvalue), ID)) +
  geom_col(aes(y = Description, fill = ONTOLOGY), alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(x = -GeneRatio, y = ID, group = 1), orientation = "y") +
  geom_point(aes(x = -GeneRatio, y = ID, colour = ONTOLOGY), size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("#A4DDD3","#F59B7B","#A8D3A0","#AC99D2")) +
  geom_text(aes(x = 0, y = Description, label = Description), 
            hjust = 0, size = 3, lineheight = 0.7) +
  facet_wrap("ONTOLOGY", scales = "free") +
  scale_fill_manual(values = c("#A4DDD3","#F59B7B","#A8D3A0","#AC99D2")) +
  labs(x = "GeneRatio(%)  and  -Log10(FDR)", y = NULL) +
  scale_x_continuous(limits = c(-max(data$GeneRatio),
                                max(-log10(data$qvalue)) + 5)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid = element_blank())

ggsave(filename="4GOKEGG.pdf",height=7,width=10)





