library(Seurat)


counts <- Read10X(data.dir = "GSM7980868_LPTC-1")
seurat <- CreateSeuratObject(counts, project="LPTC-1")



seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

library(patchwork)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)#🦋

seurat <- NormalizeData(seurat)

seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2


seurat <- ScaleData(seurat)


seurat <- RunPCA(seurat, npcs = 50)
ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))#4*6
PCHeatmap(seurat, dims = 1:10, cells = 500, balanced = TRUE, ncol = 4)#6*10


seurat <- RunTSNE(seurat, dims = 1:10)
plot1 <- TSNEPlot(seurat)
seurat <- RunUMAP(seurat, dims = 1:10)
plot2 <- UMAPPlot(seurat)
plot1 + plot2
plot2


plot2 <- FeaturePlot(seurat, c("TSC22D1","SLC34A2","MET","PDE5A","APLP2","PCSK1N","CTSB","APOC1","SDC4",'DUSP6' ,'IGSF1' ,'NUPR1' ,'NELL2' ,'GGCT' ,'DUSP5' ,'MALL' ,'PROS1'),
                     ncol=3, reduction = "umap")
ggsave("PTC.pdf", plot2, width = 10, height = 12, dpi = 300)

plot2 <- FeaturePlot(seurat, c('VCAN' ,'GNG11' ,'MAP1B' ,'ZEB2' ,'TGFB1' ,'TGFBI' ,'MMP2' ,'MMP14' ,'SERPIE1Z' ,'EB1' ,'TWIST1' ,'TWIST2'),
                     ncol=3, reduction = "umap")
ggsave("mesen.pdf", plot2, width = 9, height = 8, dpi = 300)

plot2 <- FeaturePlot(seurat, c('SORBS2' ,'SLC26A7' ,'MT1G' ,'S100A11' ,'SORD' ,'RPS2'),
                     ncol=3, reduction = "umap")
ggsave("TFC.pdf", plot2, width = 8, height = 4, dpi = 300)

plot2 <- FeaturePlot(seurat, c('EPCAM' ,'KRT8' ,'KRT18' ,'KRT7'),
                     ncol=3, reduction = "umap")
ggsave("epi.pdf", plot2, width = 8, height = 4, dpi = 300)

plot2 <- FeaturePlot(seurat, c('CD4' ,'CD8A' ,'CD14' ,'CD19','IGHG1','FGFBP2','TPSAB1','RAMP2','DCN','PDGFRB','PLP1'),
                     ncol=3, reduction = "umap")
ggsave("immunocell.pdf", plot2, width = 10, height = 10, dpi = 300)


seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 1)#[0.1,1]数字越大越精细

plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot2
ggsave("cluster.pdf", plot2, width = 6, height = 4, dpi = 300)


library(ggplot2)
ct_markers <- c('SORBS2' ,'SLC26A7' ,'MT1G' ,'S100A11' ,'SORD' ,'RPS2' ,
                'TSC22D1' ,'SLC34A2' , 'PDE5A' ,'APLP2' ,'PCSK1N' ,'CTSB' ,'APOC1' ,'SDC4' ,'DUSP6' ,'IGSF1' ,'NUPR1' ,'NELL2' ,'GGCT' ,'DUSP5' ,'MALL' ,'PROS1',
                'EPCAM' ,'KRT8' ,'KRT18' ,'KRT7' ,
                'VCAN' ,'GNG11' ,'MAP1B' ,'ZEB2' ,'TGFB1' ,'TGFBI' ,'MMP2' ,'MMP14' ,'SERPIE1Z' ,'EB1' ,'TWIST1' ,'TWIST2',
                'CD4' ,'CD8A' ,'CD14' ,'CD19','IGHG1','FGFBP2','TPSAB1','RAMP2','DCN','PDGFRB','PLP1'
                )
p=DoHeatmap(seurat, features = ct_markers) + NoLegend()
ggsave("marker.png", plot =p, width = 8, height = 6, dpi = 300)
cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p=DoHeatmap(seurat, features = top10_cl_markers$gene) + NoLegend()
ggsave("heat.png", plot =p, width = 15, height = 10, dpi = 300)

plot1 <- FeaturePlot(seurat, c('GRIN2A' ,'LRAT' ,'TYMS'), ncol = 1)
plot2 <- VlnPlot(seurat, features = c('GRIN2A' ,'LRAT' ,'TYMS'), pt.size = 0)
plot1 + plot2 + plot_layout(widths = c(1, 2))
ggsave("9.pdf", plot1 + plot2 + plot_layout(widths = c(1, 2)), width = 7, height = 10, dpi = 300)


seurat_dorsal <- subset(seurat, subset = RNA_snn_res.1 %in% c(0,1,3,4,7))
seurat_dorsal <- FindVariableFeatures(seurat_dorsal, nfeatures = 2000)
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), unlist(cc.genes))

seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c( 'TSC22D1' ,'SLC34A2' , 'PDE5A' ,'APLP2' ,'PCSK1N' ,'CTSB' ,'APOC1' ,'SDC4' ,'DUSP6' ,'IGSF1' ,'NUPR1' ,'NELL2' ,'GGCT' ,'DUSP5' ,'MALL' ,'PROS1'), ncol = 4)

library(destiny)
dm <- DiffusionMap(Embeddings(seurat_dorsal, "pca")[,1:20])
dpt <- DPT(dm)
seurat_dorsal$dpt <- rank(dpt$dpt)
FeaturePlot(seurat_dorsal, c( 'TSC22D1' ,'SLC34A2' , 'PDE5A' ,'APLP2' ,'PCSK1N' ,'CTSB' ,'APOC1' ,'SDC4' ,'DUSP6' ,'IGSF1' ,'NUPR1' ,'NELL2' ,'GGCT' ,'DUSP5' ,'MALL' ,'PROS1'), ncol=4)


if (is(seurat_dorsal[['RNA']], 'Assay5')){
  expr <- LayerData(seurat_dorsal, assay = "RNA", layer = "data")
} else{
  expr <- seurat_dorsal[['RNA']]@data
}


library(ggplot2)
plot1 <- qplot(seurat_dorsal$dpt, as.numeric(expr["VCAN",]),
               xlab="Dpt", ylab="Expression", main="VCAN") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot2 <- qplot(seurat_dorsal$dpt, as.numeric(expr["GNG11",]),
               xlab="Dpt", ylab="Expression", main="GNG11") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot3 <- qplot(seurat_dorsal$dpt, as.numeric(expr["MAP1B",]),
               xlab="Dpt", ylab="Expression", main="MAP1B") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot1 + plot2 + plot3
ggsave("mesen-time.pdf", plot1 + plot2 + plot3, width = 12, height = 7, dpi = 300)

saveRDS(seurat,file="seurat_obj_all.rds")
saveRDS(seurat_dorsal,file="seurat_obj_dorsal.rds")

