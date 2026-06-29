library(pheatmap) 
library(ComplexHeatmap)
expFile="total.normalize.txt"    
cliFile="06.clinical477.txt"
outFile="heatmap.pdf"      
var="Risk"                 


rt=read.table(expFile, sep="\t", header=T, row.names=1, check.names=F)       
Type=read.table(cliFile, sep="\t", header=T, row.names=1, check.names=F)    
gene=read.table("AI_1.3Genes.txt",sep="\t", header=T, check.names=F)
gene=gene$x
rt=rt[gene,]
rt=as.data.frame(rt)
Riskgroup=read.table("totalRisk.txt", header=T, sep="\t", check.names=F, row.names=1)  
Riskgroup=Riskgroup[rownames(Type),]
Type=cbind(Type,Riskgroup$Risk)
colnames(Type)[colnames(Type) == "Riskgroup$Risk"] <- "Risk"
Type=Type[,-grep("Tcategory2|TNM|ATA|MACIS|EORTC",colnames(Type))] 
Type=as.data.frame(Type)

Type=Type[order(Type[,var]),]   
rt=rt[,row.names(Type)]

p_num=c()
p_sci=c()
Sum=Type
exp=t(rt)
Sum=cbind(exp,Sum)
cli=Sum[,-32]
for (cli in colnames(Sum[,1:ncol(Sum)])) {
  data=Sum[c("Risk", cli)]
  colnames(data)=c("Risk", "cli")
  data=data[(data[,"cli"]!="NA"),]
  tableStat=table(data)  
  stat=chisq.test(tableStat,simulate.p.value=TRUE)   
  p <- stat[["p.value"]]
  p1 <- ifelse(p<0.001,format(p, scientific = TRUE,digits = 3), sprintf("%.03f", p))
  p2 <- ifelse(p<0.001,"<0.001", sprintf("%.03f", p)) 
  p_num[cli] <- p2
  p_sci[cli] <- p1
}

star <- ifelse(p_num<0.001,"***",ifelse(p_num<0.01,"**",ifelse(p_num<0.05,"*","")))
p_star <- paste0(p_sci," ", star)
write.table(p_star,"heatp.txt",sep="\t",quote=F,col.names = NA)

annotation_col <- data.frame(
  VitalStatus=factor(Type$VitalStatus,levels = c('Disease-specific death','Any cause death','Alive')),
  NewTumorEvent=factor(Type$NewTumorEvent),
  RAIResponse=factor(Type$RAIResponse),
  GeneFusion=factor(Type$GeneFusion),
  TERT=factor(Type$TERT),
  RAS=factor(Type$RAS),
  BRAF=factor(Type$BRAF),
  MLNCount=factor(Type$MLNCount),
  TumorDiameter=factor(Type$TumorDiameter),
  Multifocality=factor(Type$Multifocality),
  ETE=factor(Type$ETE),
  HistologicalType=factor(Type$HistologicalType),
  Mcategory=factor(Type$Mcategory),
  Ncategory = factor(Type$Ncategory),
  Tcategory = factor(Type$Tcategory),
  CombinedDisease = factor(Type$CombinedDisease),
  RadiationHistory = factor(Type$RadiationHistory, levels = c("No", "Yes", "NA")),
  Gender = factor(Type$Gender),
  Age = factor(Type$Age),
  Risk=factor(Type$Risk,levels = c("High", "Low")))
rownames(annotation_col) <- rownames(Type)
annotation_col<-annotation_col[order(annotation_col$Risk),]

col = list(Risk=c('High'="#EB7E60","Low"="#8FB4DC"),
           Age = c("<55" = "#f4e4eb", ">=55" = "#d87070"),
           Gender = c("Female" = "#d87070", "Male" = "#f4e4eb"),
           RadiationHistory=c('Yes' = "#d87070",'No'= '#f4e4eb'),
           CombinedDisease=c('Nodular hyperplasia'='#d87070','Lymphocytic thyroiditis'='#fdbfca','None'='ghostwhite','Thyroid dysfuntion'='#f4e4eb'),
           Tcategory=c('T1'="ghostwhite","T2"="#d9eed3",'T3'="#aad09d",'T4'="#62aa67"),
           Ncategory=c('N0'="#d9eed3",'N1'="#62aa67"),
           Mcategory=c('M0'="#d9eed3",'M1'="#62aa67"),
           HistologicalType=c('Classical'='ghostwhite','Follicular variant'='#aad09d','Aggressive variants'="#62aa67"),
           ETE=c('None'="ghostwhite",'Minimal '="#d9eed3",'Gross'="#62aa67"),
           Multifocality=c("Unifocal"="#d9eed3",'Multifocal'="#62aa67"),
           TumorDiameter=c('0.1~1.0'="#d9eed3","1.1~4.0"="#aad09d",'4.1~max'="#62aa67"),
           MLNCount=c('0~05'="ghostwhite","06~10"="#d9eed3",'11~20'="#aad09d",'21~max'='#62aa67'),
           BRAF=c('Wild-type'='aliceblue','Mutation'='#537cb0'),
           RAS=c('Wild-type'='aliceblue','Mutation'='#537cb0'),
           TERT=c('Wild-type'='aliceblue','Mutation'='#537cb0'),
           GeneFusion=c('Absence'='aliceblue','Presence'='#537cb0'),
           RAIResponse=c('Not received'='ghostwhite','Sensitive'='#d2d4f5','Refractory'='#6e348c'),
           NewTumorEvent=c('None'='ghostwhite','New primary tumor'='#62aa67','Distant metastasis'='#d87070','Recurrence'='#6e348c'),
           VitalStatus=c('Alive'='ghostwhite','Disease-specific death'='#6e348c','Any cause death'='#d2bcde'))

colors_na <- c("gray95")


pdf(outFile,height=8,width=15)
pheatmap(rt, annotation=Type, 
         annotation_col = annotation_col,
         annotation_colors=col,
         na_col = colors_na,  
         cluster_cols =F,    
         cluster_rows =F,
         scale = "row",
         clustering_distance_rows = "correlation",
         show_colnames=F,
         fontsize=10,
         fontsize_row=10,
         fontsize_col=8,
         cellwidth = 1.2, cellheight = 20
         )
dev.off()

