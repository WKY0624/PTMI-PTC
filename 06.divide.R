library(recipes)
library(caret)

data <- read.table("06.clinical477.txt",header=T, sep="\t", check.names=F)
data = data[,-21]
DSS <- read.table("06.time_DSS.txt",header=T, sep="\t", check.names=F)


Cli_total = merge(data,DSS, by.x = "Sample", by.y = "ID")
sam <- createDataPartition(Cli_total$Event, p = 0.6, list = F)  

trainSet <- Cli_total[sam,]
testSet <- Cli_total[-sam,]



prop.table(table(trainSet$Event))
prop.table(table(testSet$Event))
prop.table(table(Cli_total$Event))
write.table(trainSet, file = "trainSet.txt",sep="\t",quote=F,row.names=F,col.names=T)
write.table(testSet, file = "testSet.txt",sep="\t",quote=F,row.names=F,col.names=T)
write.table(Cli_total, "Cli_total.txt",,sep="\t",quote=F,row.names=F,col.names=T)


library(table1)
library(dplyr)

train <- read.table('trainSet.txt', header = T, sep = '\t', check.names = F)
test <- read.table('testSet.txt',  header = T, sep = '\t', check.names = F)
total <- read.table('Cli_total.txt', header = T, sep = '\t', check.names = F)

train$Group = c('Training cohort')  
test$Group = c('Testing cohort')    
total2 <- rbind(train,test)          

total2$DSS <-total2$DSS/365

total2$Event = ifelse(total2$Event=='1','Yes','No')

units(total2$DSS) <- "years"　

label(total2$Event)<- 'Progression events'
label(total2$TNM) <- 'AJCC 8th TNM'
label(total2$Tcategory) <- 'T category'
label(total2$Ncategory) <- 'N category'
label(total2$Mcategory) <- 'M category'

table1(~ Age + Gender + Tcategory + Ncategory + Mcategory + TNM + Event + DSS| Group, data = total2, overall = "Entire cohort")

pvalue <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    p <- t.test(y ~ g)$p.value    
  } else {
    p <- chisq.test(table(y, g))$p.value  
  }
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
}

table1(~ Age + Gender + Tcategory + Ncategory + Mcategory + TNM + Event + DSS| Group, data=total2, 
       extra.col=list(`P-value`=pvalue),  
       overall=F)
