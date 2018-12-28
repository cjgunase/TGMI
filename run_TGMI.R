#Author: Chathura J Gunasekara
#Version: 1
#data: June 20, 2017
#install.packages("doMC")
#install.packages("infotheo")
#install.packages(stats)
#install.packages(MASS)
#install.packages(iterators)
#install.packages(parallel)
#install.packages(foreach)
library("parallel")
library(stats)
library(MASS)
library(iterators)
library(parallel)
library(foreach)
library(doMC)
library(infotheo)
#Change following paramenters for input datasets
###############################################################
ncore = detectCores()

#pathway gene file name, the file must be .csv file
#sample input files are shown

                                                                                    #Pathway FILE

#Arabidopsis Lignin pathway
#pathway_file = "./input_data/Arabidopsis/pathways/lignin_pathway_data.csv"

#Arabidopsis Lignin pathway
#pathway_file = "./input_data/Arabidopsis/pathways/lignin_pathway_data.csv

#Mouse Pluripotency pathway
pathway_file = "./input_data/Mouse_pluripotency_maintenance/mouse_pathway_data.csv"

                                                                                    #TF File
#TF file name, the file must be .csv file
#Arabidopsis
#tf_file = "./input_data/Arabidopsis/Athaliana_AllTFs.csv"

#Mouse
#tf_file = "input_data/Mouse_pluripotency_maintenance/Dataset1/dataset1_135TFs.csv"
#tf_file = "input_data/Mouse_pluripotency_maintenance/Dataset2/dataset2_235TFs.csv"
tf_file = "input_data/Mouse_pluripotency_maintenance/Dataset3/dataset3_335TFs.csv"



#the format of gene expression data "samples_in_row" or "samples_in_column"
format = "samples_in_row"


#cutoff_pvalues
alpha1 = 0.05
alpha2 = 0.05

#output files
output_file1 = "./output_results/output_Ranked_TF_Frequency.csv"
output_file2 = "./output_results/output_network.csv"
output_file3 = "./output_results/output_combinatorial_TF.txt"

##############################################################



if (format == "samples_in_row") {
  pathway.data = read.csv(pathway_file,header = TRUE)
  tf.data = read.csv(tf_file,header = TRUE)
}else{
  pathway.data = read.csv(
    pathway_file,stringsAsFactors = FALSE,header = FALSE,row.names = 1
  )
  tf.data = read.csv(
    tf_file,stringsAsFactors = FALSE,header = FALSE,row.names = 1
  )
  pathway.data = t(as.matrix(pathway.data))
  tf.data = t(as.matrix(tf.data))
}

gn <-colnames(pathway.data)
tfn<-colnames(tf.data)
nsample=dim(tf.data)[1]
source("TGMI.R")
#X4CL1
#PAL1
#MYB85
#pw1<-as.numeric(pathway.data["X4CL1"]$X4CL1)
#pw2<-as.numeric(pathway.data["PAL1"]$PAL1)
#tf <-as.numeric(tf.data["MYB85"]$MYB85)
#func(pw1,pw2,tf)
results<-parallel_approach1(pathway.data,tf.data,gn,tfn)
write.csv(results,"./temp_output_all_combinations.csv",row.names = F)

###Output combinations Done###

results <- read.csv("./temp_output_all_combinations.csv",header=T)
results<-results[with(results, order(S7_div_123_pval)), ]
corr.pval<-p.adjust(results$S7_div_123_pval, method = "BH", n = length(results$S7_div_123_pval))
results$corrected.pval<-corr.pval
selected_results <-results[as.numeric(results$corrected.pval) < alpha1,]

pair_list<-rbind(data.frame(tf=selected_results$X,pw=selected_results$Y1,str=as.numeric(as.character(selected_results$S7_div_123))),
                 data.frame(tf=selected_results$X,pw=selected_results$Y2,str=as.numeric(as.character(selected_results$S7_div_123))))
#aggregate pair_list df to get average str for each pair
pair_list <- aggregate(pair_list[,3], list(pair_list$tf,pair_list$pw), mean)
colnames(pair_list) <-c("tf","pw","str")

TF_rank_MI<-data.frame(sort(table(as.character(pair_list$tf)),decreasing = T))
discardTFs<-setdiff(as.character(colnames(tf.data)),as.character(TF_rank_MI$Var1))
Freq<-rep(0,length(discardTFs))
Var1<-discardTFs
TF_rank_MI<-rbind(TF_rank_MI,data.frame(Var1,Freq))
#output 1
write.csv(TF_rank_MI,output_file1)

#output2
write.csv(pair_list,output_file2,row.names = F)

#output3
fn <- output_file3
if (file.exists(fn)) file.remove(fn)
p.vartable2TF1PW<-create_empty_table(0,4)
for (pw in as.vector(unique(pair_list$pw))){
  #pw = "CCoAOMT1"
  temp_df<-pair_list[pair_list$pw==pw,]
  if(length(unique(temp_df$tf))<2){next}
  TF2PW<-data.frame(t(combn(unique(as.character(temp_df$tf)),2)))
  pw_df<-data.frame(rep(pw,dim(TF2PW)[1]))
  colnames(pw_df)<-"pw"
  TF2PW<- cbind(pw_df,TF2PW)
  selected_2TF1PW_p.value<-c()
  for(i in 1:dim(TF2PW)[1]){
    I2<-mi3(tf.data[as.character(TF2PW[i,]$X1)],tf.data[as.character(TF2PW[i,]$X2)],
            as.vector(pathway.data[as.character(TF2PW[i,]$pw)]))
    z<-(I2-mean(randomized_3GI_2TF_1PW))/sd(randomized_3GI_2TF_1PW)
    selected_2TF1PW_p.value<-c(selected_2TF1PW_p.value,pnorm(z,lower.tail=FALSE))
    #selected_2TF1PW_p.value<- generate_permut_pval(tf.data[as.character(TF2PW[i,]$X1)],tf.data[as.character(TF2PW[i,]$X2)],
    #                     as.vector(pathway.data[as.character(TF2PW[i,]$pw)]),I2)
    
    
  }
  selected2TF1PW<-cbind(TF2PW,
                        data.frame(p.adjust(selected_2TF1PW_p.value, method = "bonferroni", n = length(selected_2TF1PW_p.value))))
  colnames(selected2TF1PW)[4]<-"corr.p.val"
  selected2TF1PW<-selected2TF1PW[with(selected2TF1PW, order(corr.p.val)), ]
  selected2TF1PW<-selected2TF1PW[1:floor(dim(selected2TF1PW)[1]*alpha2),]
  p.vartable2TF1PW<-rbind(p.vartable2TF1PW,selected2TF1PW)
  allselectedTFs<-c(as.character(selected2TF1PW$X1),as.character(selected2TF1PW$X2))
  allselectedTFs<-as.data.frame(table(allselectedTFs))
  allselectedTFs<-allselectedTFs[with(allselectedTFs,order(-Freq)),]
  catenate<-c()
  for(i in 1:dim(allselectedTFs)[1]){
    catenate<-c(catenate,paste(as.character(allselectedTFs[i,1]),"-(",allselectedTFs[i,2],")",sep = ""))
  }
                 
  
  #catenate<-unique(c(as.character(selected2TF1PW$X1),as.character(selected2TF1PW$X2)))
  catenate <- catenate[!is.na(catenate)]
  
  cat(catenate,"=>",pw,"\n",
      file="./output_results/output_combinatorial_TF.txt",
      sep=" ",append=TRUE)
  
}
#write.csv(p.vartable2TF1PW,"TFcombinationp.values.csv")
fn <- "./temp_output_all_combinations.csv"
if (file.exists(fn)) file.remove(fn)
