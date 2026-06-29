library(tidyverse)
library(ggrepel)
library(ggfun)
library(grid)
library(readxl)

df<-read.table("M1_diff_1.3MRGsall_stat.xls",row.names = 1)
colnames(df)=df[1,]
df=df[2:nrow(df),]
df$log2FoldChange <- as.numeric(df$log2FoldChange)
df$padj <- as.numeric(df$padj)
df$new_column<-rownames(df)
colnames(df)[which(colnames(df) == "new_column")] <- "ID"

p<-ggplot(data=df,
            aes(x=log2FoldChange,
                y=-log10(padj)))+
  geom_point(alpha=0.6,aes(size=-log10(padj),color=log2FoldChange))+
  scale_color_gradientn(colours=c("#6CA3D4","#7AC3DF","#FFDD8E","#F5AA61","#EB7E60"),
                        values=seq(0,1,0.2))+
  geom_point(alpha=0.7,data=df%>%
               dplyr::arrange(padj,desc(log2FoldChange))%>%
               dplyr::slice_head(n=7),
             aes(x = log2FoldChange, y = -log10(padj), size = -log10(padj)), 
             shape = 16, color = "#EB7E60", fill = "#DB6A68") + 
  geom_text_repel(data=df%>%
                    dplyr::arrange(padj,desc(log2FoldChange))%>%
                    dplyr::slice_head(n=7),
                  aes(x = log2FoldChange, y = -log10(padj), label = ID),
                  nudge_x = 0.5,
                  box.padding = 0.5,
                  nudge_y = 1,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20
  ) + 
  scale_size(range = c(2,10),
             guide = guide_legend(override.aes = list(fill = NA)))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) + 
  xlim(c(-3, 7)) + 
  geom_vline(xintercept = c(-1.5, 1.5), lty = 4, col = "black", lwd = 0.8) + 
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) + 
  xlab('log2 fold change')+
  ylab('-log10 padj')+
  theme_bw() + 
  theme(plot.title = element_text(size = 15,hjust = 0.5),
        legend.background = element_roundrect(color = '#808080',linetype = 1),
        axis.text = element_text(size = 12.5, color = "#000000"),
        axis.title = element_text(size = 15, color = "#000000")
  ) + 
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#8FB4DC")
    ), 
    xmin = -0.3, 
    xmax = -2.8,
    ymin = 96,
    ymax = 96
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Down",
      gp = grid::gpar(col = "#8FB4DC")
    ),
    xmin = -0.5, 
    xmax = -2.5,
    ymin = 98,
    ymax = 98
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
      gp = grid::gpar(lwd = 3, col = "#E1807E")
    ), 
    xmin = 3.5, 
    xmax = 5.5,
    ymin = 96,
    ymax = 96
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Up",
      gp = grid::gpar(col = "#E1807E")
    ),
    xmin = 3.5, 
    xmax = 5.5,
    ymin = 98,
    ymax = 98
  ) 
p

ggsave(filename = "1.3volcano.pdf",
       plot = p,
       height = 5,
       width = 6.5)
