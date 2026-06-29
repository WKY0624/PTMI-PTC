library(BSgenome)
library(TCGAmutations)  
library(maftools) 
library(dplyr)


tcga_available()
proj='TCGA-THCA'
laml = tcga_load(study = "THCA") 
laml = tcga_thca
gene=read.table("M1_diff_1.3MRGsall_stat.xls",header=T,sep="\t",comment.char="",check.names=F)
gene=as.character(gene$ID)

if (as.numeric(dev.cur()) != 1) graphics.off()
plotmafSummary(maf= laml, rmOutlier = TRUE,showBarcodes = FALSE,addStat = 'median', dashboard = TRUE, titvRaw = FALSE)  

oncoplot(maf = laml, top = 30, fontSize = 0.7)
colnames(laml@clinical.data)[67] = "stage"  
oncoplot(maf = laml,sortByAnnotation = TRUE, clinicalFeatures = c("stage",'gender'))


g = c('TYMS',
      'APOE',
      'ELOVL6',
      'GRIN2A',
      'GRIN2D',
      'H2BC12',
      'H2BC8',
      'H4C9',
      'LRAT',
      'MGLL',
      'NOS3',
      'SOD3',
      'TYMS',
      'UBE2C'
)   
oncoplot(maf = laml,genes = g, fontSize = 0.7)

maf = tmb(maf = laml,         
          captureSize = 80,logScale=F)
maf$sample <- substr(maf$Tumor_Sample_Barcode,1,16)
maf$sample <- gsub("-",".",maf$sample)






