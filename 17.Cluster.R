  library(limma)
  library(survival)
  library(ConsensusClusterPlus)
  
  expFile="M1_diff_1.3MRGsall_Exp.txt"   
  cliFile="06.time_DSS.txt"
  
  data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
  group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
  group=sapply(strsplit(group,""), "[", 1)
  group=gsub("2", "1", group)
  data=data[,group==0]
  data=t(data)
  data=avereps(data)
  data=log2(data+1)  
  
  rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data)) 
  
  genes=read.table("1.3train.uniCox.txt",sep="\t",check.names=F,header=T)
  genes=genes$id
  
  cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)  
  sameSample=intersect(row.names(data),row.names(cli))
  data=data[sameSample,genes]
  cli=cli[sameSample,]
  rt=cbind(cli,data)
  
  data2=t(data)
  
  sigGenes=c() 
  for(i in colnames(rt)[3:ncol(rt)]){
    cox=coxph(Surv(DSS,Event) ~ rt[,i], data = rt)  
    coxSummary=summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    if(coxP<0.05){ sigGenes=c(sigGenes,i) }
  }
  sigGenes
  data2=t(data[,sigGenes])
  
    results=ConsensusClusterPlus(
    data2,
    maxK = 9,    
    pFeature = 1,   
    title = 'PDF_1.3clusterall_log',  
    clusterAlg = "km", 
    distance = "euclidean",  
    plot = "pdf",
    writeTable = F, 
    verbose = F  
  )
  Kvec = 2:9
  x1 = 0.1; x2 = 0.9         
  PAC = rep(NA,length(Kvec)) 
  names(PAC) = paste("K=",Kvec,sep="")         
  for(i in Kvec){
    M = results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }                                 
  optK = Kvec[which.min(PAC)] 
  optK
  
  
  clusterNum = optK    
  Cluster = results[[clusterNum]][["consensusClass"]]
  Cluster = as.data.frame(Cluster)
  Cluster[,1]=paste0("C", Cluster[,1])
  ClusterOut=rbind(ID=colnames(Cluster), Cluster)
  write.table(ClusterOut, file="1.3clusterall_log.txt", sep="\t", quote=F, col.names=F)
  
