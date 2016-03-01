# stacking8. Using rank percentile ( rank/row length).  
# Based on stacking7, extract distinct molecule for each target, which yield the situation that I test the transfer learning 
# for new target and new molecule.

library(AUC)
library(data.table)
library(ggplot2)
library(gridExtra)
setwd("~/Desktop/study/SmallMolecular/stacking")
docking <- read.csv("~/Desktop/study/SmallMolecular/stacking/docking_data_dude14_20160111.txt", sep=" ",colClasses = c(rep("character",171)))
molid = docking$molid

head(docking)


########## functions
rank_percentile_na = function(x){
  na_index = which(is.na(x))
  rankPercentile = frank(x,ties.method = "average")/length(x)
  rankPercentile[na_index] = -100
  return(rankPercentile)
}

Get_feature_oriMol = function(i){
  index_target = c(1:14)
  index_feature = c(0:13) 
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  
  # target = sapply(docking$classifications,function(x) substr(x,index_target[i],index_target[i]))
  
  # cat("table for target",paste0("feature_",target_name$name)[i],":\n",table(target))
  feature = docking[,grep(paste0("X",index_feature[i],"_"),colnames(docking))]
  feature = sapply(feature,as.numeric)
  feature = apply(feature,2,function(x) x[index[[i]]])
  cat("\ndim for features",paste0("feature_",target_name$name)[i],":\n",dim(feature))
  
  feature = apply(feature,2,rank_percentile_na)
  feature = data.frame(feature)
  feature$mean_scoreType0 = apply(feature[,1:6],1,mean)  
  feature$mean_scoreType1 = apply(feature[,7:12],1,mean)  
  feature$mean_scoreTypeAll = apply(feature[,1:12],1,mean)  
  return(feature)
}

boost_cv = function(TarName){
  sink(paste0('bst_',TarName,'_8_5foldcv.txt'))
  #sink('bst_ace_3_5foldcv.txt')  
  time_start=proc.time()
  param <- list(  objective           = "binary:logistic", 
                  booster = "gbtree",
                  eval_metric = "auc",
                  eta                 = 0.04, # 0.023
                  max_depth           = 7, #changed from default of 8
                  subsample           = 0.83,  
                  colsample_bytree    = 0.77, 
                  num_parallel_tree   = 1,
                  min_child_weight    = 1,
                  gamma               = 5
                  # alpha = 0.0001, 
                  # lambda = 1
  )
  set.seed(1120)
  xgb.cv(param, eval(parse(text = paste0("dtrain_",TarName,"_v8"))),
         nrounds = 1500, 
         nfold=5,
         metrics={'error'}, 
         verbose = 1, 
         print.every.n = 10,
         maximize = TRUE,
         nthread = 12
  )
  time = time_start-proc.time()
  print(time)
  sink()
}

boost_cv_holdout_8 = function(TarName,Rounds){
  sink(paste0('bst_',TarName,'_8_5foldcv2.txt'))
  #sink('bst_ace_3_5foldcv.txt')  
  time_start=proc.time()
  param <- list(  objective           = "binary:logistic", 
                  booster = "gbtree",
                  eval_metric = "auc",
                  eta                 = 0.04, # 0.023
                  max_depth           = 7, #changed from default of 8
                  subsample           = 0.83,  
                  colsample_bytree    = 0.77, 
                  num_parallel_tree   = 1,
                  min_child_weight    = 1,
                  gamma               = 5
                  # alpha = 0.0001, 
                  # lambda = 1
  )
  set.seed(1120)
  bst_5cv = xgb.cv(param, eval(parse(text = paste0("dtrain_",TarName,"_v8"))),
                   nrounds = Rounds, 
                   nfold=5,
                   metrics={'error'}, 
                   verbose = 1, 
                   print.every.n = 10,
                   maximize = TRUE,
                   prediction = TRUE,
                   nthread = 12
  )
  time = time_start-proc.time()
  print(time)
  sink()
  bst_holdout = data.frame(Holdout_Pred = bst_5cv$pred)
  write.csv(bst_holdout, file =  paste0("bst_",TarName,"_8_5cv_holdout.csv"),row.names = F)
}


boost_train = function(TarName,Rounds){
  sink(paste0('bst_',TarName,'_8_out.txt'))
  time_start=proc.time()
  param <- list(  objective           = "binary:logistic", 
                  booster = "gbtree",
                  eval_metric = "auc",
                  eta                 = 0.04, # 0.023
                  max_depth           = 7, #changed from default of 8
                  subsample           = 0.83,  
                  colsample_bytree    = 0.77, 
                  num_parallel_tree   = 1,
                  min_child_weight    = 1,
                  gamma               = 5
                  # alpha = 0.0001, 
                  # lambda = 1
  )
  set.seed(1120)
  bst      <- xgb.train(   params              = param, 
                           data                = eval(parse(text = paste0("dtrain_",TarName,"_v8"))), 
                           nrounds             = Rounds, 
                           verbose             = 1,  #1
                           early.stop.round    = 100,
                           watchlist           = eval(parse(text = paste0("watchlist_",TarName,"_v8"))),
                           maximize            = T,
                           print.every.n = 10,
                           nthread             = 6
  )  
  
  time = time_start-proc.time()
  print(time)
  xgb.save(bst, paste0("bst_",TarName,"_8"))
  sink()
}

pred_newTarget = function(TarName){
  pred_auc = list()
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  sink(paste0("pred_auc_F",TarName,"_8.txt"))
  for( i in 1:14){
    pred_target = predict(eval(parse(text = paste0("bst_",TarName,"_8"))),data.matrix(feature_v8_oriMol[[i]]))
    pred_auc[[i]]=auc(roc(pred_target ,as.factor(target_oriMol[[i]]))) 
    print(pred_auc[[i]])
  }
  return(pred_auc)
  sink()
}

# 14 predictions from one target, total 14 targets. store in a list
pred_target = list()
for( i in 1:14){
  temp = list()
  for( j in 1:14){
    temp[[j]] = predict(eval(parse(text = paste0("bst_",target_name$name[j],"_8"))),data.matrix(feature_v8_oriMol[[i]]))
  }
  pred_target[[i]] = temp
  print(str(pred_target))
  print(i)
}

# get targets
index_target = c(1:14)
target = list()
for ( i in 1:14){
  target[[i]] = as.numeric(sapply(docking$classifications,function(x) substr(x,index_target[i],index_target[i])))
  cat("\ntable for target",paste0("feature_",target_name$name)[i],":\n",table(target[[i]]))
}


# 
original_molecule = function(i){
target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
index = grep(paste0("*_",target_name$name[i]),molid)
return(index)
}


index = list()
index[[1]] = original_molecule(1)
index[[2]] = original_molecule(2)
index[[3]] = original_molecule(3)
index[[4]] = original_molecule(4)
index[[5]] = original_molecule(5)
index[[6]] = original_molecule(6)
index[[7]] = original_molecule(7)
index[[8]] = original_molecule(8)
index[[9]] = original_molecule(9)
index[[10]] = original_molecule(10)
index[[11]] = original_molecule(11)
index[[12]] = original_molecule(12)
index[[13]] = original_molecule(13)
index[[14]] = original_molecule(14)

# generate original unqiue molecule for each target
target_oriMol = list()
target_oriMol[[1]]  = target[[1]][index[[1]]]
target_oriMol[[2]]  = target[[2]][index[[2]]]
target_oriMol[[3]]  = target[[3]][index[[3]]]
target_oriMol[[4]]  = target[[4]][index[[4]]]
target_oriMol[[5]]  = target[[5]][index[[5]]]
target_oriMol[[6]]  = target[[6]][index[[6]]]
target_oriMol[[7]]  = target[[7]][index[[7]]]
target_oriMol[[8]]  = target[[8]][index[[8]]]
target_oriMol[[9]]  = target[[9]][index[[9]]]
target_oriMol[[10]]  = target[[10]][index[[10]]]
target_oriMol[[11]]  = target[[11]][index[[11]]]
target_oriMol[[12]]  = target[[12]][index[[12]]]
target_oriMol[[13]]  = target[[13]][index[[13]]]
target_oriMol[[14]]  = target[[14]][index[[14]]]

k =1
table(target_oriMol[[k]])[2]/length(target_oriMol[[k]])
k =k+1

sink("dude14_oriMol_HitRate.txt")
k =1
table(target_oriMol[[k]])
k =k+1
sink()

# get features
feature_v8_oriMol = list()
feature_v8_oriMol[[1]] = Get_feature_oriMol(1)
feature_v8_oriMol[[2]] = Get_feature_oriMol(2)
feature_v8_oriMol[[3]] = Get_feature_oriMol(3)
feature_v8_oriMol[[4]] = Get_feature_oriMol(4)
feature_v8_oriMol[[5]] = Get_feature_oriMol(5)
feature_v8_oriMol[[6]] = Get_feature_oriMol(6)
feature_v8_oriMol[[7]] = Get_feature_oriMol(7)
feature_v8_oriMol[[8]] = Get_feature_oriMol(8)
feature_v8_oriMol[[9]] = Get_feature_oriMol(9)
feature_v8_oriMol[[10]] = Get_feature_oriMol(10)
feature_v8_oriMol[[11]] = Get_feature_oriMol(11)
feature_v8_oriMol[[12]] = Get_feature_oriMol(12)
feature_v8_oriMol[[13]] = Get_feature_oriMol(13)
feature_v8_oriMol[[14]] = Get_feature_oriMol(14)

#feature_v8_oriMol[[1]]  = apply(feature_v8[[1]],2,function(x) x[index[[1]]])


summary(feature_v8[[1]])
summary(feature_v8_oriMol[[1]])
# get xgboost dataset
dtrain_ace_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[1]]),label=as.numeric(target_oriMol[[1]]))
watchlist_ace_v8<-list(train=dtrain_ace_v8)
xgb.DMatrix.save(dtrain_ace_v8, "dtrain_ace_v8.buffer")

dtrain_adrb1_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[2]]),label=as.numeric(target_oriMol[[2]]))
watchlist_adrb1_v8<-list(train=dtrain_adrb1_v8)
xgb.DMatrix.save(dtrain_adrb1_v8, "dtrain_adrb1_v8.buffer")

dtrain_braf_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[3]]),label=as.numeric(target_oriMol[[3]]))
watchlist_braf_v8<-list(train=dtrain_braf_v8)
xgb.DMatrix.save(dtrain_braf_v8, "dtrain_braf_v8.buffer")

dtrain_cdk2_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[4]]),label=as.numeric(target_oriMol[[4]]))
watchlist_cdk2_v8<-list(train=dtrain_cdk2_v8)
xgb.DMatrix.save(dtrain_cdk2_v8, "dtrain_cdk2_v8.buffer")

dtrain_drd3_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[5]]),label=as.numeric(target_oriMol[[5]]))
watchlist_drd3_v8<-list(train=dtrain_drd3_v8)
xgb.DMatrix.save(dtrain_drd3_v8, "dtrain_drd3_v8.buffer")

dtrain_esr1_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[6]]),label=as.numeric(target_oriMol[[6]]))
watchlist_esr1_v8<-list(train=dtrain_esr1_v8)
xgb.DMatrix.save(dtrain_esr1_v8, "dtrain_esr1_v8.buffer")

dtrain_fa10_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[7]]),label=as.numeric(target_oriMol[[7]]))
watchlist_fa10_v8<-list(train=dtrain_fa10_v8)
xgb.DMatrix.save(dtrain_fa10_v8, "dtrain_fa10_v8.buffer")

dtrain_fabp4_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[8]]),label=as.numeric(target_oriMol[[8]]))
watchlist_fabp4_v8<-list(train=dtrain_fabp4_v8)
xgb.DMatrix.save(dtrain_fabp4_v8, "dtrain_fabp4_v8.buffer")

dtrain_gria2_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[9]]),label=as.numeric(target_oriMol[[9]]))
watchlist_gria2_v8<-list(train=dtrain_gria2_v8)
xgb.DMatrix.save(dtrain_gria2_v8, "dtrain_gria2_v8.buffer")

dtrain_hdac8_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[10]]),label=as.numeric(target_oriMol[[10]]))
watchlist_hdac8_v8<-list(train=dtrain_hdac8_v8)
xgb.DMatrix.save(dtrain_hdac8_v8, "dtrain_hdac8_v8.buffer")

dtrain_mmp13_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[11]]),label=as.numeric(target_oriMol[[11]]))
watchlist_mmp13_v8<-list(train=dtrain_mmp13_v8)
xgb.DMatrix.save(dtrain_mmp13_v8, "dtrain_mmp13_v8.buffer")

dtrain_pde5a_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[12]]),label=as.numeric(target_oriMol[[12]]))
watchlist_pde5a_v8<-list(train=dtrain_pde5a_v8)
xgb.DMatrix.save(dtrain_pde5a_v8, "dtrain_pde5a_v8.buffer")

dtrain_ptn1_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[13]]),label=as.numeric(target_oriMol[[13]]))
watchlist_ptn1_v8<-list(train=dtrain_ptn1_v8)
xgb.DMatrix.save(dtrain_ptn1_v8, "dtrain_ptn1_v8.buffer")

dtrain_src_v8<-xgb.DMatrix(data=data.matrix(feature_v8_oriMol[[14]]),label=as.numeric(target_oriMol[[14]]))
watchlist_src_v8<-list(train=dtrain_src_v8)
xgb.DMatrix.save(dtrain_src_v8, "dtrain_src_v8.buffer")


# load xgboost dataset
library(xgboost)
dtrain_ace_v8 <- xgb.DMatrix("dtrain_ace_v8.buffer")
watchlist_ace_v8 <- list(train=dtrain_ace_v8)

library(xgboost)
dtrain_adrb1_v8 <- xgb.DMatrix("dtrain_adrb1_v8.buffer")
watchlist_adrb1_v8 <- list(train=dtrain_adrb1_v8)

library(xgboost)
dtrain_braf_v8 <- xgb.DMatrix("dtrain_braf_v8.buffer")
watchlist_braf_v8 <- list(train=dtrain_braf_v8)

library(xgboost)
dtrain_cdk2_v8 <- xgb.DMatrix("dtrain_cdk2_v8.buffer")
watchlist_cdk2_v8 <- list(train=dtrain_cdk2_v8)

library(xgboost)
dtrain_drd3_v8 <- xgb.DMatrix("dtrain_drd3_v8.buffer")
watchlist_drd3_v8 <- list(train=dtrain_drd3_v8)

library(xgboost)
dtrain_esr1_v8 <- xgb.DMatrix("dtrain_esr1_v8.buffer")
watchlist_esr1_v8 <- list(train=dtrain_esr1_v8)

library(xgboost)
dtrain_fa10_v8 <- xgb.DMatrix("dtrain_fa10_v8.buffer")
watchlist_fa10_v8 <- list(train=dtrain_fa10_v8)

library(xgboost)
dtrain_fabp4_v8 <- xgb.DMatrix("dtrain_fabp4_v8.buffer")
watchlist_fabp4_v8 <- list(train=dtrain_fabp4_v8)

library(xgboost)
dtrain_gria2_v8 <- xgb.DMatrix("dtrain_gria2_v8.buffer")
watchlist_gria2_v8 <- list(train=dtrain_gria2_v8)

library(xgboost)
dtrain_hdac8_v8 <- xgb.DMatrix("dtrain_hdac8_v8.buffer")
watchlist_hdac8_v8 <- list(train=dtrain_hdac8_v8)

library(xgboost)
dtrain_mmp13_v8 <- xgb.DMatrix("dtrain_mmp13_v8.buffer")
watchlist_mmp13_v8 <- list(train=dtrain_mmp13_v8)

library(xgboost)
dtrain_pde5a_v8 <- xgb.DMatrix("dtrain_pde5a_v8.buffer")
watchlist_pde5a_v8 <- list(train=dtrain_pde5a_v8)

library(xgboost)
dtrain_ptn1_v8 <- xgb.DMatrix("dtrain_ptn1_v8.buffer")
watchlist_ptn1_v8 <- list(train=dtrain_ptn1_v8)

library(xgboost)
dtrain_src_v8 <- xgb.DMatrix("dtrain_src_v8.buffer")
watchlist_src_v8 <- list(train=dtrain_src_v8)


## model cv

boost_cv("ace") #[1490]	train-auc:0.997819+0.000507	train-error:0.005221+0.000401	test-auc:0.950538+0.012366	test-error:0.011784+0.001758
boost_cv("adrb1")
boost_cv("braf")
boost_cv("cdk2")
boost_cv("drd3")
boost_cv("esr1")
boost_cv("fa10")
boost_cv("fabp4")
boost_cv("gria2")
boost_cv("hdac8")
boost_cv("mmp13")
boost_cv("pde5a")
boost_cv("ptn1")
boost_cv("src")

## train model 
boost_train("ace",1350)
boost_train("adrb1",310)
boost_train("braf",1480)
boost_train("cdk2",920)
boost_train("drd3",270)
boost_train("esr1",400)
boost_train("fa10",1490)
boost_train("fabp4",1310)
boost_train("gria2",670)
boost_train("hdac8",1120)
boost_train("mmp13",1470)
boost_train("pde5a",940)
boost_train("ptn1",1320)
boost_train("src",340)



# load model
library(xgboost)
bst_ace_8 = xgb.load("bst_ace_8")
bst_adrb1_8 = xgb.load("bst_adrb1_8")
bst_braf_8 = xgb.load("bst_braf_8")
bst_cdk2_8 = xgb.load("bst_cdk2_8")
bst_drd3_8 = xgb.load("bst_drd3_8")
bst_esr1_8 = xgb.load("bst_esr1_8")
bst_fa10_8 = xgb.load("bst_fa10_8")
bst_fabp4_8 = xgb.load("bst_fabp4_8")
bst_gria2_8 = xgb.load("bst_gria2_8")
bst_hdac8_8 = xgb.load("bst_hdac8_8")
bst_mmp13_8 = xgb.load("bst_mmp13_8")
bst_pde5a_8 = xgb.load("bst_pde5a_8")
bst_ptn1_8 = xgb.load("bst_ptn1_8")
bst_src_8 = xgb.load("bst_src_8")

## predict new target

pred_auc_8 = list()
pred_auc_8[[1]] = pred_newTarget("ace")
pred_auc_8[[2]] = pred_newTarget("adrb1")
pred_auc_8[[3]] = pred_newTarget("braf")
pred_auc_8[[4]] = pred_newTarget("cdk2")
pred_auc_8[[5]] = pred_newTarget("drd3")
pred_auc_8[[6]] = pred_newTarget("esr1")
pred_auc_8[[7]] = pred_newTarget("fa10")
pred_auc_8[[8]] = pred_newTarget("fabp4")
pred_auc_8[[9]] = pred_newTarget("gria2")
pred_auc_8[[10]] = pred_newTarget("hdac8")
pred_auc_8[[11]] = pred_newTarget("mmp13")
pred_auc_8[[12]] = pred_newTarget("pde5a")
pred_auc_8[[13]] = pred_newTarget("ptn1")
pred_auc_8[[14]] = pred_newTarget("src")



##                                      # equal weight blending
#                       # mean: equal weight blending models for ace
pred_ace_all_8 = cbind(unlist(pred_target[[1]][2]),unlist(pred_target[[1]][3]),unlist(pred_target[[1]][4]),
                       unlist(pred_target[[1]][5]),unlist(pred_target[[1]][6]),unlist(pred_target[[1]][7]),
                       unlist(pred_target[[1]][8]),unlist(pred_target[[1]][9]),unlist(pred_target[[1]][10]),
                       unlist(pred_target[[1]][11]),unlist(pred_target[[1]][12]),unlist(pred_target[[1]][13]),
                       unlist(pred_target[[1]][14]))
pred_ace_mean_8 = apply(pred_ace_all_8,1,mean)
auc(roc(pred_ace_mean_8,as.factor(target_oriMol[[1]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for adrb1
pred_adrb1_all_8 = cbind(unlist(pred_target[[2]][1]),unlist(pred_target[[2]][3]),unlist(pred_target[[2]][4]),
                         unlist(pred_target[[2]][5]),unlist(pred_target[[2]][6]),unlist(pred_target[[2]][7]),
                         unlist(pred_target[[2]][8]),unlist(pred_target[[2]][9]),unlist(pred_target[[2]][10]),
                         unlist(pred_target[[2]][11]),unlist(pred_target[[2]][12]),unlist(pred_target[[2]][13]),
                         unlist(pred_target[[2]][14]))
pred_adrb1_mean_8 = apply(pred_adrb1_all_8,1,mean)
auc(roc(pred_adrb1_mean_8,as.factor(target_oriMol[[2]])))


##                                      # equal weight blending
#                       # mean: equal weight blending models for braf
pred_braf_all_8 = cbind(unlist(pred_target[[3]][1]),unlist(pred_target[[3]][2]),unlist(pred_target[[3]][4]),
                        unlist(pred_target[[3]][5]),unlist(pred_target[[3]][6]),unlist(pred_target[[3]][7]),
                        unlist(pred_target[[3]][8]),unlist(pred_target[[3]][9]),unlist(pred_target[[3]][10]),
                        unlist(pred_target[[3]][11]),unlist(pred_target[[3]][12]),unlist(pred_target[[3]][13]),
                        unlist(pred_target[[3]][14]))
pred_braf_mean_8 = apply(pred_braf_all_8,1,mean)
auc(roc(pred_braf_mean_8,as.factor(target_oriMol[[3]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for cdk2
pred_cdk2_all_8 = cbind(unlist(pred_target[[4]][1]),unlist(pred_target[[4]][2]),unlist(pred_target[[4]][3]),
                        unlist(pred_target[[4]][5]),unlist(pred_target[[4]][6]),unlist(pred_target[[4]][7]),
                        unlist(pred_target[[4]][8]),unlist(pred_target[[4]][9]),unlist(pred_target[[4]][10]),
                        unlist(pred_target[[4]][11]),unlist(pred_target[[4]][12]),unlist(pred_target[[4]][13]),
                        unlist(pred_target[[4]][14]))
pred_cdk2_mean_8 = apply(pred_cdk2_all_8,1,mean)
auc(roc(pred_cdk2_mean_8,as.factor(target_oriMol[[4]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for drd3
pred_drd3_all_8 = cbind(unlist(pred_target[[5]][1]),unlist(pred_target[[5]][2]),unlist(pred_target[[5]][3]),
                        unlist(pred_target[[5]][4]),unlist(pred_target[[5]][6]),unlist(pred_target[[5]][7]),
                        unlist(pred_target[[5]][8]),unlist(pred_target[[5]][9]),unlist(pred_target[[5]][10]),
                        unlist(pred_target[[5]][11]),unlist(pred_target[[5]][12]),unlist(pred_target[[5]][13]),
                        unlist(pred_target[[5]][14]))
pred_drd3_mean_8 = apply(pred_drd3_all_8,1,mean)
auc(roc(pred_drd3_mean_8,as.factor(target_oriMol[[5]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for esr1
pred_esr1_all_8 = cbind(unlist(pred_target[[6]][1]),unlist(pred_target[[6]][2]),unlist(pred_target[[6]][3]),
                        unlist(pred_target[[6]][4]),unlist(pred_target[[6]][5]),unlist(pred_target[[6]][7]),
                        unlist(pred_target[[6]][8]),unlist(pred_target[[6]][9]),unlist(pred_target[[6]][10]),
                        unlist(pred_target[[6]][11]),unlist(pred_target[[6]][12]),unlist(pred_target[[6]][13]),
                        unlist(pred_target[[6]][14]))
pred_esr1_mean_8 = apply(pred_esr1_all_8,1,mean)
auc(roc(pred_esr1_mean_8,as.factor(target_oriMol[[6]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for fa10
pred_fa10_all_8 = cbind(unlist(pred_target[[7]][2]),unlist(pred_target[[7]][3]),unlist(pred_target[[7]][4]),
                        unlist(pred_target[[7]][5]),unlist(pred_target[[7]][6]),unlist(pred_target[[7]][1]),
                        unlist(pred_target[[7]][8]),unlist(pred_target[[7]][9]),unlist(pred_target[[7]][10]),
                        unlist(pred_target[[7]][11]),unlist(pred_target[[7]][12]),unlist(pred_target[[7]][13]),
                        unlist(pred_target[[7]][14]))
pred_fa10_mean_8 = apply(pred_fa10_all_8,1,mean)
auc(roc(pred_fa10_mean_8,as.factor(target_oriMol[[7]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for fabp4
pred_fabp4_all_8 = cbind(unlist(pred_target[[8]][2]),unlist(pred_target[[8]][3]),unlist(pred_target[[8]][4]),
                         unlist(pred_target[[8]][5]),unlist(pred_target[[8]][6]),unlist(pred_target[[8]][7]),
                         unlist(pred_target[[8]][1]),unlist(pred_target[[8]][9]),unlist(pred_target[[8]][10]),
                         unlist(pred_target[[8]][11]),unlist(pred_target[[8]][12]),unlist(pred_target[[8]][13]),
                         unlist(pred_target[[8]][14]))
pred_fabp4_mean_8 = apply(pred_fabp4_all_8,1,mean)
auc(roc(pred_fabp4_mean_8,as.factor(target_oriMol[[8]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for gria2
pred_gria2_all_8 = cbind(unlist(pred_target[[9]][2]),unlist(pred_target[[9]][3]),unlist(pred_target[[9]][4]),
                         unlist(pred_target[[9]][5]),unlist(pred_target[[9]][6]),unlist(pred_target[[9]][7]),
                         unlist(pred_target[[9]][8]),unlist(pred_target[[9]][1]),unlist(pred_target[[9]][10]),
                         unlist(pred_target[[9]][11]),unlist(pred_target[[9]][12]),unlist(pred_target[[9]][13]),
                         unlist(pred_target[[9]][14]))
pred_gria2_mean_8 = apply(pred_gria2_all_8,1,mean)
auc(roc(pred_gria2_mean_8,as.factor(target_oriMol[[9]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for hdac8
pred_hdac8_all_8 = cbind(unlist(pred_target[[10]][2]),unlist(pred_target[[10]][3]),unlist(pred_target[[10]][4]),
                         unlist(pred_target[[10]][5]),unlist(pred_target[[10]][6]),unlist(pred_target[[10]][7]),
                         unlist(pred_target[[10]][8]),unlist(pred_target[[10]][9]),unlist(pred_target[[10]][1]),
                         unlist(pred_target[[10]][11]),unlist(pred_target[[10]][12]),unlist(pred_target[[10]][13]),
                         unlist(pred_target[[10]][14]))
pred_hdac8_mean_8 = apply(pred_hdac8_all_8,1,mean)
auc(roc(pred_hdac8_mean_8,as.factor(target_oriMol[[10]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for mmp13
pred_mmp13_all_8 = cbind(unlist(pred_target[[11]][2]),unlist(pred_target[[11]][3]),unlist(pred_target[[11]][4]),
                         unlist(pred_target[[11]][5]),unlist(pred_target[[11]][6]),unlist(pred_target[[11]][7]),
                         unlist(pred_target[[11]][8]),unlist(pred_target[[11]][9]),unlist(pred_target[[11]][10]),
                         unlist(pred_target[[11]][1]),unlist(pred_target[[11]][12]),unlist(pred_target[[11]][13]),
                         unlist(pred_target[[11]][14]))
pred_mmp13_mean_8 = apply(pred_mmp13_all_8,1,mean)
auc(roc(pred_mmp13_mean_8,as.factor(target_oriMol[[11]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for pde5a
pred_pde5a_all_8 = cbind(unlist(pred_target[[12]][2]),unlist(pred_target[[12]][3]),unlist(pred_target[[12]][4]),
                         unlist(pred_target[[12]][5]),unlist(pred_target[[12]][6]),unlist(pred_target[[12]][7]),
                         unlist(pred_target[[12]][8]),unlist(pred_target[[12]][9]),unlist(pred_target[[12]][10]),
                         unlist(pred_target[[12]][11]),unlist(pred_target[[12]][1]),unlist(pred_target[[12]][13]),
                         unlist(pred_target[[12]][14]))
pred_pde5a_mean_8 = apply(pred_pde5a_all_8,1,mean)
auc(roc(pred_pde5a_mean_8,as.factor(target_oriMol[[12]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for ptn1
pred_ptn1_all_8 = cbind(unlist(pred_target[[13]][2]),unlist(pred_target[[13]][3]),unlist(pred_target[[13]][4]),
                        unlist(pred_target[[13]][5]),unlist(pred_target[[13]][6]),unlist(pred_target[[13]][7]),
                        unlist(pred_target[[13]][8]),unlist(pred_target[[13]][9]),unlist(pred_target[[13]][10]),
                        unlist(pred_target[[13]][11]),unlist(pred_target[[13]][12]),unlist(pred_target[[13]][1]),
                        unlist(pred_target[[13]][14]))
pred_ptn1_mean_8 = apply(pred_ptn1_all_8,1,mean)
auc(roc(pred_ptn1_mean_8,as.factor(target_oriMol[[13]])))

##                                      # equal weight blending
#                       # mean: equal weight blending models for src
pred_src_all_8 = cbind(unlist(pred_target[[14]][2]),unlist(pred_target[[14]][3]),unlist(pred_target[[14]][4]),
                       unlist(pred_target[[14]][5]),unlist(pred_target[[14]][6]),unlist(pred_target[[14]][7]),
                       unlist(pred_target[[14]][8]),unlist(pred_target[[14]][9]),unlist(pred_target[[14]][10]),
                       unlist(pred_target[[14]][11]),unlist(pred_target[[14]][12]),unlist(pred_target[[14]][13]),
                       unlist(pred_target[[14]][1]))
pred_src_mean_8 = apply(pred_src_all_8,1,mean)
auc(roc(pred_src_mean_8,as.factor(target_oriMol[[14]])))

# store blending probability result
store_blending_prob = function(i){
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  write.csv(eval(parse(text = paste0("pred_",target_name$name[i],"_mean_8"))),paste0("pred_",target_name$name[i],"_mean_8.csv"), row.names = F)
  
} 
store_blending_prob(1)
store_blending_prob(2)
store_blending_prob(3)
store_blending_prob(4)
store_blending_prob(5)
store_blending_prob(6)
store_blending_prob(7)
store_blending_prob(8)
store_blending_prob(9)
store_blending_prob(10)
store_blending_prob(11)
store_blending_prob(12)
store_blending_prob(13)
store_blending_prob(14)

# roc plot
i = 1
plot(roc(eval(parse(text = paste0("pred_",target_name$name[i],"_mean_8"))) ,as.factor(target_oriMol[[i]])))
i=i+1
roc(eval(parse(text = paste0("pred_",target_name$name[i],"_mean_8"))) ,as.factor(target_oriMol[[i]]))[[3]][
  max(which(roc(eval(parse(text = paste0("pred_",target_name$name[i],"_mean_8"))) ,as.factor(target_oriMol[[i]]))[[2]]<0.01))
  ]
i=i+1

# holdout prediction
boost_cv_holdout_8("ace",1350)
boost_cv_holdout_8("adrb1",310)
boost_cv_holdout_8("braf",1480)
boost_cv_holdout_8("cdk2",920)
boost_cv_holdout_8("drd3",270)
boost_cv_holdout_8("esr1",400)
boost_cv_holdout_8("fa10",1490)
boost_cv_holdout_8("fabp4",1310)
boost_cv_holdout_8("gria2",670)
boost_cv_holdout_8("hdac8",1120)
boost_cv_holdout_8("mmp13",1470)
boost_cv_holdout_8("pde5a",940)
boost_cv_holdout_8("ptn1",1320)
boost_cv_holdout_8("src",340)


# load holdout prediction
bst_ace_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_ace_8_5cv_holdout.csv", sep="")
bst_adrb1_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_adrb1_8_5cv_holdout.csv", sep="")
bst_braf_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_braf_8_5cv_holdout.csv", sep="")
bst_cdk2_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_cdk2_8_5cv_holdout.csv", sep="")
bst_drd3_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_drd3_8_5cv_holdout.csv", sep="")
bst_esr1_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_esr1_8_5cv_holdout.csv", sep="")
bst_fa10_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_fa10_8_5cv_holdout.csv", sep="")
bst_fabp4_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_fabp4_8_5cv_holdout.csv", sep="")
bst_gria2_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_gria2_8_5cv_holdout.csv", sep="")
bst_hdac8_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_hdac8_8_5cv_holdout.csv", sep="")
bst_mmp13_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_mmp13_8_5cv_holdout.csv", sep="")
bst_pde5a_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_pde5a_8_5cv_holdout.csv", sep="")
bst_ptn1_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_ptn1_8_5cv_holdout.csv", sep="")
bst_src_8_5cv_holdout <- read.csv("~/Desktop/study/SmallMolecular/stacking/bst_src_8_5cv_holdout.csv", sep="")

# load blending prediction
pred_ace_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_ace_mean_8.csv", sep="")
pred_adrb1_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_adrb1_mean_8.csv", sep="")
pred_braf_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_braf_mean_8.csv", sep="")
pred_cdk2_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_cdk2_mean_8.csv", sep="")
pred_drd3_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_drd3_mean_8.csv", sep="")
pred_esr1_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_esr1_mean_8.csv", sep="")
pred_fa10_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_fa10_mean_8.csv", sep="")
pred_fabp4_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_fabp4_mean_8.csv", sep="")
pred_gria2_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_gria2_mean_8.csv", sep="")
pred_hdac8_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_hdac8_mean_8.csv", sep="")
pred_mmp13_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_mmp13_mean_8.csv", sep="")
pred_pde5a_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_pde5a_mean_8.csv", sep="")
pred_ptn1_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_ptn1_mean_8.csv", sep="")
pred_src_mean_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/pred_src_mean_8.csv", sep="")

# generate rank probability distribution 

get_holdoutAndTarget_8 = function(i){
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
#  paste0("bst_",target_name$name[i],"_8_5cv_hold")
 # eval(parse(text = paste0("bst_",target_name$name[i],"_8_5cv_holdout")))
  holdAndTarget = data.frame(Holdout_Pred = eval(parse(text = paste0("bst_",target_name$name[i],"_8_5cv_holdout"))),
                             target = target_oriMol[[i]]
  )
  return(holdAndTarget)
                               
}
get_blendingAndTarget_8 = function(i){
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  #  paste0("bst_",target_name$name[i],"_8_5cv_hold")
  # eval(parse(text = paste0("bst_",target_name$name[i],"_8_5cv_holdout")))
  blendingAndTarget = data.frame( Pred = eval(parse(text = paste0("pred_",target_name$name[i],"_mean_8"))),
                             target = target_oriMol[[i]]
  )
  colnames(blendingAndTarget)[1] = "Pred"
  return(blendingAndTarget)
}
df_ace_8_pred_holdout = get_holdoutAndTarget_8(1)
df_adrb1_8_pred_holdout = get_holdoutAndTarget_8(2)
df_braf_8_pred_holdout = get_holdoutAndTarget_8(3)
df_cdk2_8_pred_holdout = get_holdoutAndTarget_8(4)
df_drd3_8_pred_holdout = get_holdoutAndTarget_8(5)
df_esr1_8_pred_holdout = get_holdoutAndTarget_8(6)
df_fa10_8_pred_holdout = get_holdoutAndTarget_8(7)
df_fabp4_8_pred_holdout = get_holdoutAndTarget_8(8)
df_gria2_8_pred_holdout = get_holdoutAndTarget_8(9)
df_hdac8_8_pred_holdout = get_holdoutAndTarget_8(10)
df_mmp13_8_pred_holdout = get_holdoutAndTarget_8(11)
df_pde5a_8_pred_holdout = get_holdoutAndTarget_8(12)
df_ptn1_8_pred_holdout = get_holdoutAndTarget_8(13)
df_src_8_pred_holdout = get_holdoutAndTarget_8(14)

df_ace_8_pred_blending = get_blendingAndTarget_8(1)
df_adrb1_8_pred_blending = get_blendingAndTarget_8(2)
df_braf_8_pred_blending = get_blendingAndTarget_8(3)
df_cdk2_8_pred_blending = get_blendingAndTarget_8(4)
df_drd3_8_pred_blending = get_blendingAndTarget_8(5)
df_esr1_8_pred_blending = get_blendingAndTarget_8(6)
df_fa10_8_pred_blending = get_blendingAndTarget_8(7)
df_fabp4_8_pred_blending = get_blendingAndTarget_8(8)
df_gria2_8_pred_blending = get_blendingAndTarget_8(9)
df_hdac8_8_pred_blending = get_blendingAndTarget_8(10)
df_mmp13_8_pred_blending = get_blendingAndTarget_8(11)
df_pde5a_8_pred_blending = get_blendingAndTarget_8(12)
df_ptn1_8_pred_blending = get_blendingAndTarget_8(13)
df_src_8_pred_blending = get_blendingAndTarget_8(14)

hist_plot_cv_8 = function(i){
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  png(paste0("./plot/hist_holdout_", target_name$name[i],"_8.png"),units = "in", width=11, height=8.5, res=300)
  par(mfrow=c(2,1))
  hist(eval(parse(text = paste0("df_",target_name$name[i],"_8_pred_holdout")))$Holdout_Pred[which(eval(parse(text = paste0("df_",target_name$name[i],"_8_pred_holdout")))$target==1)],col = 2, ylim = c(1,60),breaks = seq(0, 1, by =0.01),
       main = paste0("Histogram of Active molecule for ",target_name$name[i]), xlab = "Predicted probability from cv Holdout")
  hist(eval(parse(text = paste0("df_",target_name$name[i],"_8_pred_holdout")))$Holdout_Pred[which(eval(parse(text = paste0("df_",target_name$name[i],"_8_pred_holdout")))$target==0)],col = 3, ylim = c(1,60),breaks = seq(0, 1, by =0.01),
       main = paste0("Histogram of Decoy molecule for ",target_name$name[i]), xlab = "Predicted probability from cv Holdout")
    dev.off()
}

hist_plot_blending_8 = function(i){
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  png(paste0("./plot/hist_blending_", target_name$name[i],"_8.png"),units = "in", width=11, height=8.5, res=300)
  par(mfrow=c(2,1))
  hist(eval(parse(text = paste0("df_",target_name$name[i],"_8_pred_blending")))$Pred[which(eval(parse(text = paste0("df_",target_name$name[i],"_8_pred_blending")))$target==1)],col = 2, ylim = c(1,60),breaks = seq(0, 1, by =0.01),
       main = paste0("Histogram of Active molecule for ",target_name$name[i]), xlab = "Predicted probability from blending")
  hist(eval(parse(text = paste0("df_",target_name$name[i],"_8_pred_blending")))$Pred[which(eval(parse(text = paste0("df_",target_name$name[i],"_8_pred_blending")))$target==0)],col = 3, ylim = c(1,60),breaks = seq(0, 1, by =0.01),
       main = paste0("Histogram of Decoy molecule for ",target_name$name[i]), xlab = "Predicted probability from cv blending")
  dev.off()
}

hist_plot_cv_8(1)
hist_plot_cv_8(2)
hist_plot_cv_8(3)
hist_plot_cv_8(4)
hist_plot_cv_8(5)
hist_plot_cv_8(6)
hist_plot_cv_8(7)
hist_plot_cv_8(8)
hist_plot_cv_8(9)
hist_plot_cv_8(10)
hist_plot_cv_8(11)
hist_plot_cv_8(12)
hist_plot_cv_8(13)
hist_plot_cv_8(14)

hist_plot_blending_8(1)
hist_plot_blending_8(2)
hist_plot_blending_8(3)
hist_plot_blending_8(4)
hist_plot_blending_8(5)
hist_plot_blending_8(6)
hist_plot_blending_8(7)
hist_plot_blending_8(8)
hist_plot_blending_8(9)
hist_plot_blending_8(10)
hist_plot_blending_8(11)
hist_plot_blending_8(12)
hist_plot_blending_8(13)
hist_plot_blending_8(14)

# write function to calculate enrichment.
get_enrichmentTable_cv_8 = function(k){
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  random_hit = length(which(target_oriMol[[k]]==1))/length(target_oriMol[[k]])
  holdout_with_target = data.frame(holdout_pred = eval(parse(text = paste0("bst_",target_name$name[k],"_8_5cv_holdout$Holdout_Pred"))),
                                   label = target_oriMol[[k]])
  holdout_with_target = holdout_with_target[order(holdout_with_target$holdout_pred,decreasing = T),]
  nfold = numeric(4999)
  j=1
  for( i in seq(1,dim(holdout_with_target)[1],dim(holdout_with_target)[1]*0.0002)[-1]){
    
    nfold[j]=(length(which(holdout_with_target$label[1:i] == 1))/i)/random_hit
    j = j+1
  }
  enrichment_table = data.frame(folds = nfold, scale = seq(0,1,0.0002)[-c(1,5001)], TarName = target_name$name[k])
  write.csv(enrichment_table, file = paste0("enrichment_",target_name$name[k],"_cv_8.csv"),row.names = F)
  return(enrichment_table)
}

get_enrichmentTable_blending_8 = function(k){
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  random_hit = length(which(target_oriMol[[k]]==1))/length(target_oriMol[[k]])
  holdout_with_target = data.frame(holdout_pred = eval(parse(text = paste0("pred_",target_name$name[k],"_mean_8"))),
                                   label = target_oriMol[[k]])
  
  holdout_with_target = holdout_with_target[order(holdout_with_target$holdout_pred,decreasing = T),]
  nfold = numeric(4999)
  j=1
  for( i in seq(1,dim(holdout_with_target)[1],dim(holdout_with_target)[1]*0.0002)[-1]){
    
    nfold[j]=(length(which(holdout_with_target$label[1:i] == 1))/i)/random_hit
    j = j+1
  }
  enrichment_table = data.frame(folds = nfold, scale = seq(0,1,0.0002)[-c(1,5001)], TarName = target_name$name[k])
  write.csv(enrichment_table, file = paste0("enrichment_",target_name$name[k],"_blending_8.csv"),row.names = F)
  return(enrichment_table)
}


# generate enrichment data file
enrichment_ace_cv_8 = get_enrichmentTable_cv_8(1)
enrichment_adrb1_cv_8 = get_enrichmentTable_cv_8(2)
enrichment_braf_cv_8 = get_enrichmentTable_cv_8(3)
enrichment_cdk2_cv_8 = get_enrichmentTable_cv_8(4)
enrichment_drd3_cv_8 = get_enrichmentTable_cv_8(5)
enrichment_esr1_cv_8 = get_enrichmentTable_cv_8(6)
enrichment_fa10_cv_8 = get_enrichmentTable_cv_8(7)
enrichment_fabp4_cv_8 = get_enrichmentTable_cv_8(8)
enrichment_gria2_cv_8 = get_enrichmentTable_cv_8(9)
enrichment_hdac8_cv_8 = get_enrichmentTable_cv_8(10)
enrichment_mmp13_cv_8 = get_enrichmentTable_cv_8(11)
enrichment_pde5a_cv_8 = get_enrichmentTable_cv_8(12)
enrichment_ptn1_cv_8 = get_enrichmentTable_cv_8(13)
enrichment_src_cv_8 = get_enrichmentTable_cv_8(14)

enrichment_ace_blending_8 = get_enrichmentTable_blending_8(1)
enrichment_adrb1_blending_8 = get_enrichmentTable_blending_8(2)
enrichment_braf_blending_8 = get_enrichmentTable_blending_8(3)
enrichment_cdk2_blending_8 = get_enrichmentTable_blending_8(4)
enrichment_drd3_blending_8 = get_enrichmentTable_blending_8(5)
enrichment_esr1_blending_8 = get_enrichmentTable_blending_8(6)
enrichment_fa10_blending_8 = get_enrichmentTable_blending_8(7)
enrichment_fabp4_blending_8 = get_enrichmentTable_blending_8(8)
enrichment_gria2_blending_8 = get_enrichmentTable_blending_8(9)
enrichment_hdac8_blending_8 = get_enrichmentTable_blending_8(10)
enrichment_mmp13_blending_8 = get_enrichmentTable_blending_8(11)
enrichment_pde5a_blending_8 = get_enrichmentTable_blending_8(12)
enrichment_ptn1_blending_8 = get_enrichmentTable_blending_8(13)
enrichment_src_blending_8 = get_enrichmentTable_blending_8(14)


# load enrichment 
enrichment_ace_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_ace_cv_8.csv")
enrichment_adrb1_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_adrb1_cv_8.csv")
enrichment_braf_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_braf_cv_8.csv")
enrichment_cdk2_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_cdk2_cv_8.csv")
enrichment_drd3_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_drd3_cv_8.csv")
enrichment_esr1_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_esr1_cv_8.csv")
enrichment_fa10_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_fa10_cv_8.csv")
enrichment_fabp4_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_fabp4_cv_8.csv")
enrichment_gria2_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_gria2_cv_8.csv")
enrichment_hdac8_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_hdac8_cv_8.csv")
enrichment_mmp13_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_mmp13_cv_8.csv")
enrichment_pde5a_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_pde5a_cv_8.csv")
enrichment_ptn1_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_ptn1_cv_8.csv")
enrichment_src_cv_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_src_cv_8.csv")
#
enrichment_ace_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_ace_blending_8.csv")
enrichment_adrb1_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_adrb1_blending_8.csv")
enrichment_braf_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_braf_blending_8.csv")
enrichment_cdk2_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_cdk2_blending_8.csv")
enrichment_drd3_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_drd3_blending_8.csv")
enrichment_esr1_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_esr1_blending_8.csv")
enrichment_fa10_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_fa10_blending_8.csv")
enrichment_fabp4_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_fabp4_blending_8.csv")
enrichment_gria2_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_gria2_blending_8.csv")
enrichment_hdac8_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_hdac8_blending_8.csv")
enrichment_mmp13_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_mmp13_blending_8.csv")
enrichment_pde5a_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_pde5a_blending_8.csv")
enrichment_ptn1_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_ptn1_blending_8.csv")
enrichment_src_blending_8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/enrichment_src_blending_8.csv")


# plot function
enrichment_plot_cv_8=function(i){
  stacking8_cv_auc=c(0.950708,0.961903,0.923995,0.935525,0.865923,0.954746,0.963578,0.885355,0.927422,0.985011,0.959252,0.941147,0.944218,0.862725)
  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  
  k<-ggplot(eval(parse(text = paste0("enrichment_",target_name$name[i],"_cv_8"))), aes(scale, folds))+
    geom_point(color="blue",alpha=1/10)+
    xlim(c(0,0.05))+ylim(c(0,1000))+labs(title=paste0("Version 8-CV Enrichment plot for ",target_name$name[i],"\nCrossValidation AUC =",stacking8_cv_auc[i]),x = "Threshold")+
    geom_hline(yintercept = 100)
  #print(k)
  return(k)
}


enrichment_plot_blending_8=function(i){
  stacking8_blending_auc=c(0.8879094,0.9291406,0.889576,0.9111473,0.8232989,0.863151,0.9003936,0.7694326,0.8738852,0.9311819,0.8929501,0.8980644,0.8541414,0.7449019)

  target_name = data.frame(name = c("ace","adrb1","braf","cdk2","drd3","esr1","fa10","fabp4","gria2","hdac8","mmp13","pde5a","ptn1","src"))
  k<-ggplot(eval(parse(text = paste0("enrichment_",target_name$name[i],"_blending_8"))), aes(scale, folds))+
    geom_point(color="blue",alpha=1/10)+
    xlim(c(0,0.05))+ylim(c(0,1000))+labs(title=paste0("Version 8-Blending Enrichment plot for ",target_name$name[i],"\nAUC =",stacking8_blending_auc[i]),x = "Threshold")+
    geom_hline(yintercept = 100)
  #print(k)
  return(k)
}

enrichment_plot_cv_8_1=enrichment_plot_cv_8(1)
enrichment_plot_cv_8_2=enrichment_plot_cv_8(2)
enrichment_plot_cv_8_3=enrichment_plot_cv_8(3)
enrichment_plot_cv_8_4=enrichment_plot_cv_8(4)
enrichment_plot_cv_8_5=enrichment_plot_cv_8(5)
enrichment_plot_cv_8_6=enrichment_plot_cv_8(6)
enrichment_plot_cv_8_7=enrichment_plot_cv_8(7)
enrichment_plot_cv_8_8=enrichment_plot_cv_8(8)
enrichment_plot_cv_8_9=enrichment_plot_cv_8(9)
enrichment_plot_cv_8_10=enrichment_plot_cv_8(10)
enrichment_plot_cv_8_11=enrichment_plot_cv_8(11)
enrichment_plot_cv_8_12=enrichment_plot_cv_8(12)
enrichment_plot_cv_8_13=enrichment_plot_cv_8(13)
enrichment_plot_cv_8_14=enrichment_plot_cv_8(14)

enrichment_plot_blending_8_1=enrichment_plot_blending_8(1)
enrichment_plot_blending_8_2=enrichment_plot_blending_8(2)
enrichment_plot_blending_8_3=enrichment_plot_blending_8(3)
enrichment_plot_blending_8_4=enrichment_plot_blending_8(4)
enrichment_plot_blending_8_5=enrichment_plot_blending_8(5)
enrichment_plot_blending_8_6=enrichment_plot_blending_8(6)
enrichment_plot_blending_8_7=enrichment_plot_blending_8(7)
enrichment_plot_blending_8_8=enrichment_plot_blending_8(8)
enrichment_plot_blending_8_9=enrichment_plot_blending_8(9)
enrichment_plot_blending_8_10=enrichment_plot_blending_8(10)
enrichment_plot_blending_8_11=enrichment_plot_blending_8(11)
enrichment_plot_blending_8_12=enrichment_plot_blending_8(12)
enrichment_plot_blending_8_13=enrichment_plot_blending_8(13)
enrichment_plot_blending_8_14=enrichment_plot_blending_8(14)


grid.arrange(enrichment_plot_cv_8_1,enrichment_plot_cv_8_2,enrichment_plot_cv_8_3,
             enrichment_plot_cv_8_4,enrichment_plot_cv_8_5,enrichment_plot_cv_8_6,
             enrichment_plot_cv_8_7,enrichment_plot_cv_8_8,enrichment_plot_cv_8_9,
             enrichment_plot_cv_8_10,enrichment_plot_cv_8_11,enrichment_plot_cv_8_12,
             enrichment_plot_cv_8_13,enrichment_plot_cv_8_14,ncol=3)

grid.arrange(enrichment_plot_blending_8_1,enrichment_plot_blending_8_2,enrichment_plot_blending_8_3,
             enrichment_plot_blending_8_4,enrichment_plot_blending_8_5,enrichment_plot_blending_8_6,
             enrichment_plot_blending_8_7,enrichment_plot_blending_8_8,enrichment_plot_blending_8_9,
             enrichment_plot_blending_8_10,enrichment_plot_blending_8_11,enrichment_plot_blending_8_12,
             enrichment_plot_blending_8_13,enrichment_plot_blending_8_14,ncol=3)


# enrichments of 14 targets together
stacking8_cv_auc=c(0.950708,0.961903,0.923995,0.935525,0.865923,0.954746,0.963578,0.885355,0.927422,0.985011,0.959252,0.941147,0.944218,0.862725)

enrichment_all_cv_8 = rbind(enrichment_ace_cv_8,enrichment_adrb1_cv_8,enrichment_braf_cv_8,
                            enrichment_cdk2_cv_8,enrichment_drd3_cv_8,enrichment_esr1_cv_8,
                            enrichment_fa10_cv_8,enrichment_fabp4_cv_8,enrichment_gria2_cv_8,
                            enrichment_hdac8_cv_8,enrichment_mmp13_cv_8,enrichment_pde5a_cv_8,
                            enrichment_ptn1_cv_8,enrichment_src_cv_8)

enrichment_plot_all_cv_8 = ggplot(enrichment_all_cv_8, aes(scale, folds,color = as.factor(TarName)))+
  geom_point(alpha=1)+
  xlim(c(0,0.05))+ylim(c(0,1000))+
  geom_hline(yintercept = 100)+
  labs(title = paste0("Version 8-cv Enrichment plot all.\nOverall AUC mean = ",mean(stacking8_cv_auc)))
print(enrichment_plot_all_cv_8)

# different ylim
enrichment_plot_all_cv_8_2 = ggplot(enrichment_all_cv_8, aes(scale, folds,color = as.factor(TarName)))+
  geom_point(alpha=1)+
  xlim(c(0,0.05))+ylim(c(0,100))+
  geom_hline(yintercept = 10)+
  labs(title = paste0("Version 8-cv Enrichment plot all.\nOverall AUC mean = ",mean(stacking8_cv_auc)))
print(enrichment_plot_all_cv_8_2)

png(paste0("enrichment_plot_all_cv_8_2.png"), units = "in", width=11, height=8.5, res=300)
dev.off()
#
stacking8_blending_auc=c(0.8879094,0.9291406,0.889576,0.9111473,0.8232989,0.863151,0.9003936,0.7694326,0.8738852,0.9311819,0.8929501,0.8980644,0.8541414,0.7449019)

enrichment_all_blending_8 = rbind(enrichment_ace_blending_8,enrichment_adrb1_blending_8,enrichment_braf_blending_8,
                                  enrichment_cdk2_blending_8,enrichment_drd3_blending_8,enrichment_esr1_blending_8,
                                  enrichment_fa10_blending_8,enrichment_fabp4_blending_8,enrichment_gria2_blending_8,
                                  enrichment_hdac8_blending_8,enrichment_mmp13_blending_8,enrichment_pde5a_blending_8,
                                  enrichment_ptn1_blending_8,enrichment_src_blending_8)

enrichment_plot_all_blending_8 = ggplot(enrichment_all_blending_8, aes(scale, folds,color = as.factor(TarName)))+
  geom_point(alpha=1)+
  xlim(c(0,0.05))+ylim(c(0,1000))+
  geom_hline(yintercept = 100)+
  labs(title = paste0("Version 8-Blending Enrichment plot all.\nOverall AUC mean = ",mean(stacking8_blending_auc)))
print(enrichment_plot_all_blending_8)

# different ylim
enrichment_plot_all_blending_8_2 = ggplot(enrichment_all_blending_8, aes(scale, folds,color = as.factor(TarName)))+
  geom_point(alpha=1)+
  xlim(c(0,0.05))+ylim(c(0,100))+
  geom_hline(yintercept = 10)+
  labs(title = paste0("Version 8-Blending Enrichment plot all.\nOverall AUC mean = ",mean(stacking8_blending_auc)))
print(enrichment_plot_all_blending_8_2)

png(paste0("enrichment_plot_all_blending_8_2.png"), units = "in", width=11, height=8.5, res=300)
dev.off()





#
a=c(0.8546438	,0.7437953,	0.8846257,	0.8341925,	0.8084892,	0.7614526,	0.7432026,	0.7887209,	0.8224245,	0.8954837,	0.8439154,	0.8157085,	0.7230949)
b = c(         
               0.0153521,
               0.01505845,
               0.01674675,
               0.01391224,
               0.01819823,
               0.02611868,
               0.01680973, 
               0.01317765,
               0.01601055,
               0.01517322,
               0.01425552, 
               0.0176319, 
               0.01498085) 
plot(a,b)
          

blend_hit_stack8 <- read.csv("~/Desktop/study/SmallMolecular/stacking/blend_hit_stack8.csv", sep=";")
head(blend_hit_stack8)

plot(blend_hit_stack8$ace[-1],blend_hit_stack8$hitRate[-1])

ggplot(blend_hit_stack8,aes(ace,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(adrb1,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(braf,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(cdk2,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(drd3,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(esr1,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(fa10,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(fabp4,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(gria2,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(hdac8,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(mmp13,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(pde5a,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(ptn1,hitRate))+geom_point()
ggplot(blend_hit_stack8,aes(,hitRate))+geom_point()