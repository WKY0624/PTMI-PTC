library(Rtsne)
library(scatterplot3d)
library(ggplot2)

bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F,row.names=1)
  data=rt[c(3:(ncol(rt)-2))]
  Risk=rt[,"Risk"]  
    data.pca=prcomp(data, scale. = TRUE)
  pcaPredict=predict(data.pca)
  group = levels(factor(Risk))
  col=mycol[match(Risk,group)]
  PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2], PC3 = pcaPredict[,3],
                   Risk=Risk)	
  pdf(file=pcaFile, height = 4.5, width = 4.5)
  p <- scatterplot3d(PCA[,c(1,2,3)], grid = T, box=T, angle = 45, pch = 19, color = alpha(col,0.5),   #train:c(7,3,5)
                     col.axis="gray30",  
                     xlab = "PC1", ylab = "PC2", zlab = "PC3"
  )
  print(p)
  dev.off()  
  
  tsneOut=Rtsne(data, dims=3, perplexity=10, verbose=F, max_iter=500, check_duplicates=F)
  tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2], tSNE3 = tsneOut$Y[,3], Risk=Risk)	
  pdf(file=tsneFile, height=4.5, width=4.5)      
  p <- scatterplot3d(tsne[,c(2,3,1)], grid = T, box=T, angle = 45, pch = 19, color = alpha(col,0.5),  #train:c(3,1,2)
                     col.axis="gray30",   
                     xlab = "tSNE1", ylab = "tSNE2", zlab = "tSNE3"
  )
  print(p)
  dev.off()
  
  
}


bioPCA(inputFile="trainRisk.txt", pcaFile="PCA.train.pdf", tsneFile="t-SNE.train.pdf")
bioPCA(inputFile="testRisk.txt", pcaFile="PCA.test.pdf", tsneFile="t-SNE.test.pdf")
bioPCA(inputFile="totalRisk.txt", pcaFile="PCA.total.pdf", tsneFile="t-SNE.total.pdf")
