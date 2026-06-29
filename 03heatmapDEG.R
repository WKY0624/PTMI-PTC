library(pheatmap)
library(dplyr)
gene=read.table("01.newCounts_477+59.txt",header=T,sep="\t",comment.char="",row.names =1)
gene50=read.table("M1_diff_1.3MRGsall_stat.xls",header=T,sep="\t",comment.char="",row.names =1)
top50_genes <- gene50 %>% arrange(padj) %>% head(50) %>% rownames()
top50_genes
heatmap_data <- gene[top50_genes, ]
p=gene50[top50_genes, ]
heatmap_data=cbind(heatmap_data,p)

group_list=read.table("01.Group.txt")
group_list=factor(group_list$condition,levels=c("Normal","Tumor"))

annotation_col<-data.frame(Group=group_list)
rownames(annotation_col)<-colnames(gene)

my_breaks <- seq(-5, 5, length.out = 101)

heatmap_plot <- pheatmap(heatmap_data,
                         annotation_col=annotation_col,
                         color = colorRampPalette(c('#537cb0',"#8FB4DC",'ghostwhite','#EB7E60','#d87070'))(100),
                         show_rownames = TRUE,
                         show_colnames = FALSE,
                         breaks = my_breaks,
                         cluster_rows = F,
                         cluster_cols = F,
                         fontsize_row = 6,
                         fontsize = 5,
                         scale = "row",
                         main = "Top 50 Differentially Expressed Genes")

pdf("Heatmap_Top50_DEGs.pdf", width = 4, height = 4)
print(heatmap_plot)
dev.off()
