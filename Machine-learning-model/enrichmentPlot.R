# Based on model_building
# Generate enrichment plot for holdout predictions and test predictions.

library(AUC)
library(data.table)
library(ggplot2)
library(gridExtra)
setwd("~/Desktop/study/SmallMolecular/stacking")
docking <- read.csv("~/Desktop/study/SmallMolecular/stacking/docking_data_dude14_20160111.txt", sep=" ",colClasses = c(rep("character",171)))
molid = docking$molid


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
