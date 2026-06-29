library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggsci)
library(survival)
library(randomForestSRC)
library(plsRcox)
library(superpc)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(ggbreak)
library(tidyr)
library(ggbreak)
library(edgeR)
library(limma)
library(survival)
library(survminer)
library(stringi)
library(ggplot2)
library(ggpubr)
library(beepr)
library(pheatmap)
library(data.table)
library(ggsignif)
library(RColorBrewer)
library(future.apply)
library(gplots)
library(DESeq2)
library(ggrepel)
library(Rcpp)
library(survivalsvm)
library(dplyr)
library(rms)
library(pec)
library(ggDCA)
library(glmnet)
library(foreign)
library(regplot)
library(randomForestSRC)
library(timeROC)
library(tidyr)
library(tibble)
library(caret)
library(gbm)
library(tidyverse)
library(obliqueRSF)
library(remotes)
library(aorsf)
library(xgboost)
library(party)
library(partykit)

uniGenes <-"1.3train.uniSigExp.txt"
unidata <- read.table(uniGenes, header = T, row.names = 1, sep = '\t',check.names = F)
genedata <- t(unidata[,-grep('Event',colnames(unidata))])

con <- unidata[unidata$Event == '0', ]
treat <- unidata[unidata$Event == '1', ]

conNum = nrow(con)
treatNum = nrow(treat)

type = c(rep("con",conNum),rep("treat",treatNum))

merge2 <- rbind(con,treat)
merge2 <- merge2[,-grep('Event',colnames(merge2))]
merge2 <- t(merge2)

out=rbind(id=paste0(colnames(merge2),"_",type), merge2)
write.table(out, file="1.3data.txt", sep="\t", quote=F, col.names=F)

data = read.table('1.3data.txt',header = T, row.names = 1, sep = '\t',check.names =F)
data = t(data)
data=data[,2:ncol(data)]
group=gsub("(.*)\\_(.*)", "\\2", rownames(data))



train<-read.table(file="1.3train_expTime_DSS.txt", header = T,sep = "\t", quote = "", check.names = F,row.names = 1)
test<-read.table(file="1.3test_expTime_DSS.txt",header = T,sep = "\t", quote = "", check.names = F,row.names = 1)

unidata=unidata[,3:ncol(unidata)]
trainset<-train[,colnames(unidata)]
trainset<-c(train[,1:2],trainset)
trainset<-as.data.frame(trainset)
row.names(trainset)=row.names(unidata)
test<-test[,colnames(trainset)]
train=trainset

str(train)
train$DSS=as.numeric(train$DSS)
train$Event=as.numeric(train$Event)

test$DSS=as.numeric(test$DSS)
test$Event=as.numeric(test$Event)

trainlist=list(Train=train,Test=test)

result <- data.frame()

#### 1-1.RSF####
set.seed(seed)
fit <- rfsrc(Surv(DSS,Event)~.,data = trainlist$Train,
             nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(DSS,Event)~.,data = trainlist$Train,
             ntree = best,nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result,cc)



#### 1-2.RSF + Enet####
vi <- data.frame(imp=vimp.rfsrc(fit)$importance)
vi$imp <- (vi$imp-min(vi$imp))/(max(vi$imp)-min(vi$imp))
vi$ID <- rownames(vi)

pdf("rsf_highgene(rsf).pdf")
ggplot(vi,aes(imp,reorder(ID,imp)))+
  geom_bar(stat = 'identity',fill='#FF9933',color='black',width=0.7)+
  geom_vline(xintercept = 0.01,color='grey50',linetype=2)+
  labs(x='Relative importance by Random Forest',y=NULL)+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.x = element_text(size = 11,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.title = element_text(size=13,color='black'),
        legend.text = element_text(size=12,color='black'),
        legend.title = element_text(size=13,color='black'))+
  scale_y_discrete(expand = c(0.03,0.03))+
  scale_x_continuous(expand = c(0.01,0.01))
dev.off()
rid <- rownames(vi)[vi$imp>0.01]
train2 <- train[,c('DSS','Event',rid)]
trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})

x1 <- as.matrix(train2[,rid])
x2 <- as.matrix(Surv(train2$DSS,train2$Event))

for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}
set.seed(seed)
modelexp=as.matrix(trainlist2$Train[,c(3:ncol(trainlist2$Train))])
modelstat=Surv(trainlist2$Train$DSS,trainlist2$Train$Event)
Enetmodel <- glmnet(modelexp,modelstat,family = 'cox',nfolds=10)
Enetmodel_cv<-cv.glmnet(modelexp,modelstat,family = 'cox',nfolds=10)
model_fit<-glmnet(modelexp,modelstat,family = 'cox',nfolds=10,keep=T,lambda = Enetmodel_cv$lambda.min)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model_fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=model_fit$lambda.min)))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + Enet','[lambda=',round(Enetmodel_cv$lambda.min,3),']')
result <- rbind(result,cc)

set.seed(seed)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + StepCox','[',direction,']')
  result <- rbind(result,cc)
}




#### 1-4.RSF + CoxBoost ####
set.seed(seed)
pen <- optimCoxBoostPenalty(train2[,'DSS'],train2[,'Event'],as.matrix(train2[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train2[,'DSS'],train2[,'Event'],as.matrix(train2[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(train2[,'DSS'],train2[,'Event'],as.matrix(train2[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBoost')
result <- rbind(result,cc)

#### 1-5.RSF + plsRcox ####
set.seed(seed)
model_exp=data.frame(train2[,-c(1:2)])
model_time=train2$DSS
model_stat=train2$Event
model<-plsRcox(model_exp,time = model_time,event = model_stat,nt=10)
cv.model<-cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = T)
model<-plsRcox(model_exp,
               time = model_time,
               event = model_stat,
               nt=cv.model$lambda.min5,
               alpha.pvals.expli = 0.05,
               sparse = T,
               pvals.expli = T)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model,type="lp",newdata=x[,-c(1,2)])))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + plsRcox')
result <- rbind(result,cc)


data <- list(x=t(train2[,-c(1,2)]),y=train2$DSS,censoring.status=train2$Event,featurenames=colnames(train2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)

rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$DSS,censoring.status=w$Event,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + SuperPC')
result <- rbind(result,cc)

#### 1-7.RSF + GBM ####
set.seed(seed)
fit <- gbm(formula = Surv(DSS,Event)~.,
           data = train2,
           distribution = 'coxph',
           n.minobsinnode = 10,
           n.cores = 1,
           shrinkage = 0.005,
           interaction.depth = 2,
           cv.folds = 5)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(DSS,Event)~.,
           data = train2,
           distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,
           n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + GBM')
result <- rbind(result,cc)

#### 1-8.RSF + survivalsvm ####
set.seed(seed)
fit = survivalsvm(Surv(DSS,Event)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + survival-SVM')
result <- rbind(result,cc)

#### 1-9.RSF + Ridge ####
set.seed(seed)
modelexp=as.matrix(train2[,c(3:ncol(train2))])
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  model <- glmnet(modelexp,train2$Event,family = 'binomial',alpha = alpha,nfolds=10)
  model_cv<-cv.glmnet(modelexp,train2$Event,family = 'binomial',alpha =alpha,nfolds=10)
  fit<-glmnet(modelexp,train2$Event,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + Ridge','[α=',alpha,']')
  result <- rbind(result,cc)
}

#### 1-10.RSF + obliqueRSF ####
set.seed(seed)
model<-orsf(data = train2,n_tree = 9,formula = Surv(DSS,Event)~.)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + obliqueRSF')
result <- rbind(result,cc)


Genes<-model[["importance"]]
Genes=as.data.frame(Genes)
Genes=rownames(Genes)
write.table(Genes, file="AI_1.3Genes.txt",sep="\t",quote=F,row.names = F)
trainrisk=rs$Train
testrisk=rs$Test
totalrisk=rbind(trainrisk,testrisk)
genes_with_importance=model$importance
write.table(genes_with_importance, file="genecoef.txt", sep="\t", quote=F, row.names=F)
write.table(trainrisk,"trainrisk.txt",sep="\t",quote=F,col.names = NA)
write.table(testrisk,"testrisk.txt",sep="\t",quote=F,col.names = NA)
write.table(totalrisk,"totalrisk.txt",sep="\t",quote=F,col.names = NA)

GenesExp=genedata[Genes,,drop=F]
write.table(GenesExp, file="AI_1.3GenesExp.xls",sep="\t",quote=F)

#### 1-11.RSF + xgboost ####
set.seed(seed)
model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$DSS)
object<-list(bojective="surivival:cox",
             booster="gbtree",
             eval_metric="cox-nloglik",
             eta=0.01,
             max_depth=3,
             subsample=1,
             colsample_bytree=1,
             gamma=0.5)
model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + xgboost')
result <- rbind(result,cc)



#### 2-1.Enet ####
modelexp=as.matrix(train[,c(3:ncol(train))])
modelstat=Surv(train$DSS,train$Event)
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  model <- glmnet(modelexp,modelstat,family = 'cox',alpha = alpha,nfolds=10)
  model_cv<-cv.glmnet(modelexp,modelstat,family = 'cox',alpha = alpha,nfolds=10)
  fit <-glmnet(modelexp,modelstat,family = 'cox',alpha =alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
  rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}

#### 2-2.Lasso + RSF####
set.seed(seed)
fit = cv.glmnet(modelexp, modelstat,family = "cox")
coef.min = coef(fit, s = "lambda.min") 
rid <- coef.min@Dimnames[[1]]

train2 <- train[,c('DSS','Event',rid)]
trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})

set.seed(seed)
fit <- rfsrc(Surv(DSS,Event)~.,data = train2,
             ntree = 1000,nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(DSS,Event)~.,data = train2,
             ntree = best,nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'Lasso + RSF'
result <- rbind(result,cc)

#### 2-3.Lasso + StepCox ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

#### 2-4.Lasso + CoxBoost ####
set.seed(seed)
modelpen<-optimCoxBoostPenalty(time = train2$DSS,
                               status = train2$Event,
                               as.matrix(train2[,-c(1:2)]),
                               trace = T,
                               parallel = T)
cvmodel<-cv.CoxBoost(time = train2$DSS,
                     status = train2$Event,
                     as.matrix(train2[,-c(1:2)]),
                     maxstepno = 100,
                     K = 3,
                     type = "verweij",
                     penalty=modelpen$penalty)
fit<-CoxBoost(time = train2$DSS,
              status = train2$Event,
              as.matrix(train2[,-c(1:2)]),
              stepno = cvmodel$optimal.step,
              penalty = modelpen$penalty)

rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result,cc)

#### 2-5.Lasso + plsRcox ####
set.seed(seed)
model_exp=data.frame(train2[,-c(1:2)])
model_time=train2$DSS
model_stat=train2$Event
model<-plsRcox(model_exp,time = model_time,event = model_stat,nt=10)
cv.model<-cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = F)
cv.plsRcox.res=cv.plsRcox(list(x=model_exp,time=model_time,status=model_stat),nt=5,verbose = F)
fit <- plsRcox(model_exp,
               time = model_time,
               event = model_stat,
               nt=cv.model$lambda.min5,
               alpha.pvals.expli = 0.05,
               sparse = T,
               pvals.expli = T)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + plsRcox')
result <- rbind(result,cc)

#### 2-6.Lasso + superpc ####
data <- list(x=t(train2[,-c(1,2)]),y=train2$DSS,censoring.status=train2$Event,featurenames=colnames(train2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$DSS,censoring.status=w$Event,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + SuperPC')
result <- rbind(result,cc)

#### 2-7.Lasso + GBM ####
set.seed(seed)
fit <- gbm(formula = Surv(DSS,Event)~.,data = train2,distribution = 'coxph',
           n.trees = 1000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(DSS,Event)~.,data = train2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + GBM')
result <- rbind(result,cc)

#### 2-8.Lasso + survivalsvm ####
set.seed(seed)
fit = survivalsvm(Surv(DSS,Event)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + survival-SVM')
result <- rbind(result,cc)

#### 2-9.Lasso + Ridge ###
set.seed(seed)
modelexp=as.matrix(train2[,c(3:ncol(train2))])
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  model <- glmnet(modelexp,train2$Event,family = 'binomial',alpha = alpha,nfolds=10)
  model_cv<-cv.glmnet(modelexp,train2$Event,family = 'binomial',alpha =alpha,nfolds=10)
  fit<-glmnet(modelexp,train2$Event,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + Ridge','[α=',alpha,']')
  result <- rbind(result,cc)
}


#### 2-10.Lasso + obliqueRSF ####
set.seed(seed)
model<-orsf(data = train2,formula = Surv(DSS,Event)~.)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + obliqueRSF')
result <- rbind(result,cc)


 #### 2-11.Lasso + xgboost ####
set.seed(seed)
model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$DSS)
object<-list(bojective="surivival:cox",
             booster="gbtree",
             eval_metric="cox-nloglik",
             eta=0.01,
             max_depth=3,
             subsample=1,
             colsample_bytree=1,
             gamma=0.5)
model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + xgboost')
result <- rbind(result,cc)


#### 3-1.StepCox ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']')
  result <- rbind(result,cc)
}

#### 3-2.StepCox + RSF ####
for (direction in c("both","backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  
  set.seed(seed)
  fit <- rfsrc(Surv(DSS,Event)~.,data = train2,
               ntree = 1000,nodesize = rf_nodesize,
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  best <- which.min(fit$err.rate)
  set.seed(seed)
  fit <- rfsrc(Surv(DSS,Event)~.,data = train2,
               ntree = best,nodesize = rf_nodesize,
               splitrule = 'logrank',
               importance = T,
               proximity = T,
               forest = T,
               seed = seed)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + RSF')
  result <- rbind(result,cc)
}

#### 3-3.StepCox + Enet ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  x1 <- as.matrix(train2[,rid])
  x2 <- as.matrix(Surv(train2$DSS,train2$Event))
  
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
    rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
    rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
    rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
    rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
    rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox','[',direction,']',' + Enet','[α=',alpha,']')
    result <- rbind(result,cc)
  }
}
#### 3-4.StepCox + CoxBoost ####
for (direction in c("both", "backward", "forward") ){
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  
  set.seed(seed)
  pen <- optimCoxBoostPenalty(train2[,'DSS'],train2[,'Event'],as.matrix(train2[,-c(1,2)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(train2[,'DSS'],train2[,'Event'],as.matrix(train2[,-c(1,2)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit <- CoxBoost(train2[,'DSS'],train2[,'Event'],as.matrix(train2[,-c(1,2)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + CoxBoost')
  result <- rbind(result,cc)
}

#### 3-5.StepCox + plsRcox ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  cv.plsRcox.res=cv.plsRcox(list(x=train2[,rid],time=train2$DSS,status=train2$Event),nt=10,nfold = 10,verbose = F)
  fit <- plsRcox(train2[,rid],time=train2$DSS,event=train2$Event,nt=as.numeric(cv.plsRcox.res[5]))
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + plsRcox')
  result <- rbind(result,cc)
}

#### 3-6.StepCox + superpc ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  
  data <- list(x=t(train2[,-c(1,2)]),y=train2$DSS,censoring.status=train2$Event,featurenames=colnames(train2)[-c(1,2)])
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                       n.fold = 5,
                       n.components=3,
                       min.features=1,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  rs <- lapply(trainlist2,function(w){
    test <- list(x=t(w[,-c(1,2)]),y=w$DSS,censoring.status=w$Event,featurenames=colnames(w)[-c(1,2)])
    ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
    rr <- as.numeric(ff$v.pred)
    rr2 <- cbind(w[,1:2],RS=rr)
    return(rr2)
  })
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + SuperPC')
  result <- rbind(result,cc)
}

#### 3-7.StepCox + gbm ####😩
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  
  set.seed(seed)
  fit <- gbm(formula = Surv(DSS,Event)~.,data = train2,distribution = 'coxph',
             n.trees = 1000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 5,n.cores = 1)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(DSS,Event)~.,data = train2,distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 5,n.cores = 1)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + GBM')
  result <- rbind(result,cc)
}


#### 3-8.StepCox + survival-SVM ####
for (direction in c("both", "backward", "forward")) {
  #direction='both'
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  
  fit = survivalsvm(Surv(DSS,Event)~., data= train2, gamma.mu = 1)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']',' + survival-SVM')
  result <- rbind(result,cc)
}

#### 3-9.StepCox + Ridge ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  set.seed(seed)
  modelexp=as.matrix(train2[,c(3:ncol(train2))])
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    model <- glmnet(modelexp,train2$Event,family = 'binomial',alpha = alpha,nfolds=10)
    model_cv<-cv.glmnet(modelexp,train2$Event,family = 'binomial',alpha =alpha,nfolds=10)
    fit<-glmnet(modelexp,train2$Event,family = 'binomial',alpha = alpha,nfolds=10,keep=T,lambda = model_cv$lambda.min)
    rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="response",newx=as.matrix(x[,-c(1,2)]))))})
    rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
    rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
    rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
    rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
    cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
      rownames_to_column('ID')
    cc$Model <- paste0('StepCox','[',direction,']','+ Ridge','[α=',alpha,']')
    result <- rbind(result,cc)
  }
}

#### 3-10.StepCox + obliqueRSF ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
  model<-orsf(data = train2,n_tree = 100,formula = Surv(DSS,Event)~.)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, new_data=x,pred_type = "risk")[,1]))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']', '+ obliqueRSF')
  result <- rbind(result,cc)
}


#### 3-11.StepCox + xgboost ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train),direction = direction)
  rid <- names(coef(fit))
  train2 <- train[,c('DSS','Event',rid)]
  trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})
    model_mat<-xgb.DMatrix(data = as.matrix(train2[,-c(1:2)]),label=train2$DSS)
  object<-list(bojective="surivival:cox",
               booster="gbtree",
               eval_metric="cox-nloglik",
               eta=0.01,
               max_depth=3,
               subsample=1,
               colsample_bytree=1,
               gamma=0.5)
  model<-xgb.train(params=object,data = model_mat,nrounds = 100,watchlist = list(val2=model_mat),early_stopping_rounds = 10)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(model, newdata=as.matrix(x[,-c(1:2)]))))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox','[',direction,']', '+ xgboost')
  result <- rbind(result,cc)
}



#### 4-1.CoxBoost ####
set.seed(seed)
pen <- optimCoxBoostPenalty(train[,'DSS'],train[,'Event'],as.matrix(train[,-c(1,2)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(train[,'DSS'],train[,'Event'],as.matrix(train[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(train[,'DSS'],train[,'Event'],as.matrix(train[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result,cc)

#### 4-2.CoxBoost + Enet ####

rid <- names(coef(fit)[which(coef(fit)!=0)])
train2 <- train[,c('DSS','Event',rid)]
trainlist2 <- lapply(trainlist,function(x){x[,c('DSS','Event',rid)]})

x1 <- as.matrix(train2[,rid])
x2 <- as.matrix(Surv(train2$DSS,train2$Event))

for (alpha in seq(0.1,1,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + Enet','[α=',alpha,']')
  result <- rbind(result,cc)
}


#### 4-3.CoxBoost + stepcox ####

for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(DSS,Event)~.,train2),direction = direction)
  rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
  rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
  rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
  cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + StepCox','[',direction,']')
  result <- rbind(result,cc)
}

#### 4-4.CoxBoost + RSF ####

set.seed(seed)
fit <- rfsrc(Surv(DSS,Event)~.,data = train2,
             ntree = 100,nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
best <- which.min(fit$err.rate)
set.seed(seed)
fit <- rfsrc(Surv(DSS,Event)~.,data = train2,
             ntree = best,nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- 'CoxBoost + RSF'
result <- rbind(result,cc)

#### 4-5.CoxBoost + plsRcox ####
cv.plsRcox.res=cv.plsRcox(list(x=train2[,rid],time=train2$DSS,status=train2$Event),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(train2[,rid],time=train2$DSS,event=train2$Event,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + plsRcox')
result <- rbind(result,cc)

#### 4-6.CoxBoost + superpc ####
data <- list(x=t(train2[,-c(1,2)]),y=train2$DSS,censoring.status=train2$Event,featurenames=colnames(train2)[-c(1,2)])
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
cv.fit <- superpc.cv(fit,data,n.threshold = 20, 
                     n.fold = 5,
                     n.components=3,
                     min.features=1,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(trainlist2,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$DSS,censoring.status=w$Event,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + SuperPC')
result <- rbind(result,cc)

#### 4-7.CoxBoost + GBM ####

set.seed(seed)
fit <- gbm(formula = Surv(DSS,Event)~.,data = train2,distribution = 'coxph',
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 3,n.cores = 1)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(DSS,Event)~.,data = train2,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 3,n.cores = 1)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + GBM')
result <- rbind(result,cc)

#### 4-8.CoxBoost + survivalsvm ####

fit = survivalsvm(Surv(DSS,Event)~., data= train2, gamma.mu = 2)
rs <- lapply(trainlist2,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + survival-SVM')
result <- rbind(result,cc)

#### 5.plsRcox####

cv.plsRcox.res=cv.plsRcox(list(x=train[,-c(1,2)],time=train$DSS,status=train$Event),nt=10,nfold = 10,verbose = F)
fit <- plsRcox(train[,-c(1,2)],time=train$DSS,event=train$Event,nt=as.numeric(cv.plsRcox.res[5]))
rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result,cc)




#### 6.superpc ####

data <- list(x=t(train[,-c(1,2)]),y=train$DSS,censoring.status=train$Event,featurenames=colnames(train)[-c(1,2)])
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) 
cv.fit <- superpc.cv(fit,data,n.threshold = 20,
                     n.fold = 10,
                     n.components=3,
                     min.features=1, 
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs <- lapply(trainlist,function(w){
  test <- list(x=t(w[,-c(1,2)]),y=w$DSS,censoring.status=w$Event,featurenames=colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit,data,test,threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('SuperPC')
result <- rbind(result,cc)

#### 7.GBM ####

set.seed(seed)
fit <- gbm(formula = Surv(DSS,Event)~.,data = train,distribution = 'coxph',
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(DSS,Event)~.,data = train,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 5,n.cores = 1)

rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result,cc)

#### 8.survivalsvm ####

fit = survivalsvm(Surv(DSS,Event)~., data= train, gamma.mu = 2)
rs <- lapply(trainlist,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="Inf",]
rs[["Train"]]=rs[["Train"]][rs[["Train"]]$RS !="-Inf",]
rs[["Test"]]=rs[["Test"]][rs[["Test"]]$RS!="-Inf",]
cc <- data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(DSS,Event)~RS,x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('survival-SVM')


result <- rbind(result,cc)
result$Cindex=round(result$Cindex,3)
result2 <- result 
result2=setDT(result2)  
result2_train=result2[result2$ID=="Train",]
result2_train=as.data.frame(result2_train)
result2_test=result2[result2$ID=="Test",]
result2_test=as.data.frame(result2_test)
Cindexnums=data.frame(Model=result2_train$Model,Train=result2_train$Cindex,Test=result2_test$Cindex)
Cindexnums[,-1] <- apply(Cindexnums[,-1], 2, as.numeric)
Cindexnums$All <- apply(Cindexnums[,2:3], 1, mean)
Cindexnums <- Cindexnums[order(Cindexnums$Test, decreasing = T),]
write.table(Cindexnums,"out_Cindex.txt", col.names = T, row.names = F, sep = "\t", quote = F)
nums <- Cindexnums[, 2:3]%>%as.matrix()
rownames(nums)=Cindexnums$Model
Cindex_mat=nums
avg_Cindex <- apply(Cindex_mat, 1, mean)     
avg_Cindex <- sort(avg_Cindex, decreasing = T)    
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]      
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) 
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                          gp = gpar(fill = "#8FB4DC", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 9, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)
names(CohortCol) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)

cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha, 
              top_annotation = col_ha,
              col = c('#8FB4DC', "#FFFFFF", '#EB7E60'), 
              rect_gp = gpar(col = "black", lwd = 1), 
              cluster_columns = FALSE, cluster_rows = FALSE, 
              show_column_names = FALSE, 
              show_row_names = TRUE,
              row_names_side = "left",
              width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
              height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { 
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf(file.path( "Cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 5, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
invisible(dev.off())




