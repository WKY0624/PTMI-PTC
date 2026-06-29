library(shapviz)
library(h2o)
library(modeldata)          
library(tidymodels)
data =read.table("Cli_total.txt",header=T,sep="\t",comment.char="",check.names=F,row.names = 1)
traindata=read.table("trainSet.txt",header=T,sep="\t",comment.char="",check.names=F,row.names = 1)
testdata=read.table("testSet.txt",header=T,sep="\t",comment.char="",check.names=F,row.names = 1)
risk=read.table("totalRisk.txt",header=T,sep="\t",comment.char="",check.names=F,row.names = 1)
data=data[,grep("Age|Gender|Tcategory|Ncategory|Mcategory|Event",colnames(data))]
data=data[,-3]
data=data[,-6]
traindata=traindata[,colnames(data)]
testdata=testdata[,colnames(data)]

risk=risk[rownames(data),]
data=cbind(risk$Risk,data)
colnames(data)[colnames(data) == "risk$Risk"] <- "Risk"

data$Age=ifelse(data$Age=='<55',0,1)
data$Gender=ifelse(data$Gender=='Male',0,1)
data$Tcategory=ifelse(data$Tcategory=='T1',1,
                       ifelse(data$Tcategory=='T2',2,
                              ifelse(data$Tcategory=='T3',3,4)))
data$Ncategory=ifelse(data$Ncategory=='N0',0,1)
data$Mcategory=ifelse(data$Mcategory=='M0',0,1)

data$Risk=ifelse(data$Risk=='High',1,0)
traindata=data[rownames(traindata),]
testdata=data[rownames(testdata),]

traindata$Event=factor(traindata$Event)
testdata$Event=factor(testdata$Event)
data$Event=factor(data$Event)

h2o.init()

traindata_env <- recipe(Event~.,traindata) %>%  
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%             
  bake(new_data=NULL) %>%            
  as.h2o() 

testdata_env <- recipe(Event~.,testdata) %>%    
  step_dummy(all_nominal_predictors()) %>%            
  prep() %>%             
  bake(new_data=NULL) %>%            
  as.h2o() 

data_env <- recipe(Event~.,data) %>%    
  step_dummy(all_nominal_predictors()) %>%            
  prep() %>%             
  bake(new_data=NULL) %>%            
  as.h2o() 

traindata_env$Age=h2o.asfactor(traindata_env$Age)
traindata_env$Gender=h2o.asfactor(traindata_env$Gender)
traindata_env$Tcategory=h2o.asfactor(traindata_env$Tcategory)
traindata_env$Ncategory=h2o.asfactor(traindata_env$Ncategory)
traindata_env$Mcategory=h2o.asfactor(traindata_env$Mcategory)
traindata_env$Risk=h2o.asfactor(traindata_env$Risk)

 testdata_env$Age=h2o.asfactor( testdata_env$Age)
 testdata_env$Gender=h2o.asfactor( testdata_env$Gender)
 testdata_env$Tcategory=h2o.asfactor( testdata_env$Tcategory)
 testdata_env$Ncategory=h2o.asfactor( testdata_env$Ncategory)
 testdata_env$Mcategory=h2o.asfactor( testdata_env$Mcategory)
 testdata_env$Risk=h2o.asfactor( testdata_env$Risk)

 data_env$Age=h2o.asfactor( data_env$Age)
 data_env$Gender=h2o.asfactor( data_env$Gender)
 data_env$Tcategory=h2o.asfactor( data_env$Tcategory)
 data_env$Ncategory=h2o.asfactor( data_env$Ncategory)
 data_env$Mcategory=h2o.asfactor( data_env$Mcategory)
  data_env$Risk=h2o.asfactor( data_env$Risk)



fit = h2o.gbm(colnames(traindata[,-4]), "Event", training_frame = traindata_env)
shap = shapviz(fit, X_pred = traindata)
sv_force(shap,
         fill_colors = c("#EB7E60","#8FB4DC"),
         row_id = 1L,
         show_annotation = TRUE,
         bar_label_size = 3)

sv_dependence(shap, 
              "Risk")
sv_importance(shap,kind = "beeswarm")+ theme_light()

fit = h2o.randomForest( x=colnames(traindata[-6]),y = "Event",                             
                        training_frame = traindata_env, nfolds = 5                           
                        )
shap = shapviz(fit, X_pred = traindata)
sv_force(shap, fill_colors = c("#EB7E60","#8FB4DC"),
         row_id = 3L,
         show_annotation = TRUE,
         bar_label_size = 3)
sv_importance(shap,kind = "beeswarm")+ theme_light()


multi_models_fit <- h2o.automl(y="Event",            
                               training_frame =traindata_env,                
                               max_models = 10
                               )
first_5_model <- h2o.get_leaderboard(multi_models_fit)
best <- h2o.get_best_model(multi_models_fit)
shap = shapviz(best, X_pred = data)

final_values <- rowSums(shap$S) + shap$baseline
shap$final_values <- final_values
shap.risk=shap$S
shap.risk=cbind(shap.risk,final_values)
write.table(shap.risk,"shap.risk.txt",sep="\t",quote=F,col.names = NA)

perf <- h2o.performance(best)
perf
h2o.auc(perf)
plot(perf,type="roc")

h2o.varimp_plot(best)
h2o.permutation_importance_plot(best,data_env)

sv_force(shap,
         fill_colors = c("#EB7E60","#8FB4DC"),
         row_id = 1L,
         show_annotation = TRUE,
         bar_label_size = 3)


sv_dependence(shap,c('Age','Gender','Tcategory','Ncategory','Mcategory'))


contributions <- h2o.predict_contributions(best, data_env)
contributions
h2o.shap_summary_plot(best, data_env)
h2o.shap_explain_row_plot(best, data_env, row_index = 1L)

coef=h2o.coef_norm(best, traindata_env) 
h2o.std_coef_plot(best) 


h2o.pd_plot(best,data_env,"Risk")
h2o.pd_multi_plot(multi_models_fit@leaderboard,data_env, "Risk")
h2o.ice_plot(best,data_env,show_pdp = TRUE,"Risk")
h2o.learning_curve_plot(best)


fit = h2o.gbm(colnames(testdata[,-7]), "Event", training_frame = testdata_env)
shap = shapviz(fit, X_pred = testdata)
sv_force(shap,
         fill_colors = c("#EB7E60","#8FB4DC"),
         row_id = 3L,
         show_annotation = TRUE,
         bar_label_size = 3)

sv_dependence(shap, 'Risk')


fit = h2o.gbm(colnames(data[,-7]), "Event", training_frame = data_env)
shap = shapviz(fit, X_pred = data)
sv_force(shap,
         fill_colors = c("#EB7E60","#8FB4DC"),
         row_id = 1L,
         show_annotation = TRUE,
         bar_label_size = 3)

sv_dependence(shap, 'Risk')
