
ccleExp <- read.table('32.OmicsExpressionProteinCodingGenesTPMLogp1(240223).csv',header = T, row.names = 1, sep = ',', check.names = F)
thyroid <- read.table("32.cell lines in Thyroid.csv",header = T, row.names = 1, sep = ',', check.names = F, quote = "")
same <- intersect(rownames(ccleExp),rownames(thyroid))
ccleTHCA <- ccleExp[same,,drop = F]
thyroid <- thyroid[same,,drop = F]
colnames(ccleTHCA) = gsub("\\(.*?\\)","",colnames(ccleTHCA))   
colnames(ccleTHCA) = as.character(colnames(ccleTHCA))
ccleTHCA = as.data.frame(ccleTHCA)

merge2 = cbind(thyroid,ccleTHCA)
merge2 <- cbind(id=row.names(merge2),merge2)
write.table(merge2, 'CCLE_THCA_TPM(240223).txt',col.names = T, row.names = F, quote = F, sep = '\t')

ccleTHCA <- read.table("CCLE_THCA_TPM(240223).txt",header = T,sep = '\t', check.names = F)

filterGene <- colnames(ccleTHCA) %in% c("Cell Line","Tumor Type","Primary Disease",
                                        'ATAD2','ELOVL6','GRIN2A','GRIN2D','H4C9','LPCAT1','NOS3','LRAT','MGLL','SOD3','TYMS','WIF1','UBE2C')
                                        

merge <- ccleTHCA[, filterGene]
merge$Abbr <- ifelse(merge$`Primary Disease`=='Well-Differentiated Thyroid Cancer','Well-DTC',
                     ifelse(merge$`Primary Disease`=='Anaplastic Thyroid Cancer','ATC',
                            ifelse(merge$`Primary Disease`=='Medullary Thyroid Cancer','MTC','Poorly-DTC')))
merge <- merge[order(merge$Abbr, decreasing =F),]


data <- data.frame(x=merge$`Cell Line`,y=merge$TYMS, tumortype=merge$`Tumor Type`, Abbr=merge$Abbr)
data$tumortype = as.factor(data$tumortype)
data$x=as.factor(data$x)

data <- data[order(data$y, decreasing = F),]   
data$x <- factor(data$x,levels = data$x)


ggplot(data, aes(x=x, y=y, color=data$Abbr)) +
  geom_segment(aes(x=x, xend=x, y=0, yend=y, linetype=factor(data$tumortype,levels = c("Primary","Metastatic"))), 
               color=ifelse(data$tumortype %in% c("Primary"), "black", "black"),
               
  ) +
  geom_point(
    size=2.5,
    alpha=0.9, pch=19) + 
  scale_color_manual(values = c("#EB7E60",'#62aa67','#d2d4f5',"#8FB4DC"), name="Primary disease") +
  scale_linetype_manual(values = c(1,2), name ='Tumor type')+
  coord_flip() +
  theme_test(base_size = 10, base_line_size = 0.4, base_rect_size = 0.5) +
  theme(axis.text = element_text(size = 8, color = 'black'),
        axis.title.y = element_text(vjust = 3),
        axis.title.x =element_text(vjust = -2),
        plot.margin = margin(10,10,10,10)) +
  xlab("Cell line") +
  ylab("TYMS expression") +
  theme(legend.position = 'right',
        legend.key.size = unit(10,'pt'),
        legend.title = element_text(size = 8, lineheight = 4),
        legend.text = element_text(size = 8, lineheight = 5),) +
  guides(linetype=guide_legend(order = 0, keyheight=0.8),
         color=guide_legend(order = 1, keyheight=0.8, override.aes = list(size=1.5))) 


ggsave('TYMS.pdf', height = 3, width = 3.7)


ref <- read.table("OmicsProfiles.csv", header = T, sep = ",", check.names = F)
ref <- ref[,c('ProfileID','ModelID')]
rownames(ref) <- ref[,1]


count <- read.table("OmicsExpressionGenesExpectedCountProfile.csv", header = T, sep = ",", check.names = F)
rownames(count) <- count[,1]
colnames(count)[1] <- "ProfileID"   

intersect <- intersect(rownames(count),rownames(ref))
count2 <- count[intersect,,drop=F]
ref2 <- ref[intersect,,drop=F]

data <- merge(ref2, count2, by.x = 'ProfileID', by.y = "ProfileID")  
data2<- subset(data,select=-c(ProfileID)) 
write.table(data2,'CCLE_counts.txt',sep="\t",quote=F, col.names=T, row.names = F)

ccleExp <- read.table('CCLE_count.txt',header = T, sep = '\t', quote = "",check.names = F)
ccleExp_mt <- as.matrix(ccleExp)
rownames(ccleExp_mt)=ccleExp_mt[,1]    
exp=ccleExp_mt[,2:ncol(ccleExp_mt)]    
dimnames=list(rownames(exp),colnames(exp))    
ccleExp_mt2=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
same <- intersect(rownames(ccleExp_mt2),rownames(thyroid))
ccleTHCA <- ccleExp_mt2[same,,drop = F]
thyroid <- thyroid[same,,drop = F]
colnames(ccleTHCA) = gsub("\\(.*?\\)","",colnames(ccleTHCA))    
colnames(ccleTHCA) = as.character(colnames(ccleTHCA))
ccleTHCA = as.data.frame(ccleTHCA)