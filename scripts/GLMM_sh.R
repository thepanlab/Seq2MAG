require("tidyverse")
require("mgcv")
require("MASS")
library("dplyr")
library("glmmTMB")
library("lme4")
library("bbmle")

args = commandArgs(trailingOnly=TRUE)

# Load meta_data
meta_data = read.csv("meta_data_shuffle.csv", head = T)
rownames(meta_data) = meta_data$Run

# Choose the k-th shuffle, k = 0 means initial result
k = as.numeric(args[1])
if(k > 0){
  meta_data = meta_data[,-which(names(meta_data)==c('pre_after_IA'))]
  colnames(meta_data)[which(colnames(meta_data)==c(paste('shuffle_IA', k, sep='_')))] = c('pre_after_IA')
}

# Read abundance data
data = read.table(args[2], head=T, sep="\t", stringsAsFactors=F)
data = t(data) #transpose data
colnames(data) = data[1,] #colnames of data are values of the first row
data = data[-1,] #remove the first row
data = as.data.frame(data) #make it a data frame

#Fixed effects and random effects
mixed_effects = c('age_at_collection', 'delivery_simple','cc','breastfeeding', 'solid_food', 'pre_after_IA', 'dbgap_maskid')
fixed_effects = mixed_effects[-which(mixed_effects == c('dbgap_maskid'))]

#result settings
all_result = c()

# Model results for all clusters

for(cluster in 1:ncol(data)){
  #Check process
  if(cluster%%10==0){
    print("The computation is coming to cluster ")
    print(cluster)
  }
  #Abundance for each cluster
  abundance = data[,cluster]
  abundance = as.data.frame(abundance)
  abundance[,1] = as.character(abundance[,1])
  abundance[,1] = as.numeric(abundance[,1])
  abundance[,1] = round(abundance[,1])
  rownames(abundance) = rownames(data)
  
  #data containing abundance
  data_glmm = meta_data[,mixed_effects]
  data_glmm[,'si'] = meta_data[,'num_wgs_reads']
  data_glmm = merge(data_glmm, abundance, by = 0)
  data_glmm = na.omit(data_glmm)
  average = mean(data_glmm$si)
  total_reads = sum(data_glmm$si)
  data_glmm$si = data_glmm$si/average
  
  #glmm formula
  glmm_formula = paste("abundance", paste('offset(log(si))', paste(fixed_effects, collapse="+"), '(1|dbgap_maskid)', sep='+'), sep="~")
  glmm_formula = as.formula(glmm_formula)
  data_glmm[,"abundance"] = as.integer(data_glmm[,"abundance"])
  
  #calculation and result summary
  glmm_model = glmmTMB(glmm_formula, data = data_glmm, family = nbinom2)
  summary = summary(glmm_model)
  #the sum of read counts
  result = as.data.frame(total_reads)
  #the number of subjects
  No_subjects = summary$ngrps$cond
  result = cbind(result, No_subjects)
  #the number of observations
  No_observations = summary$nobs
  result = cbind(result, No_observations)
  #coefficients
  coefficients = summary$coefficients$cond
  #set up colnames of coefficients
  col = c()
  for(parameter in colnames(coefficients)){
    for(feature in rownames(coefficients)){
      col = append(col, paste(feature, parameter, sep='_'))
    }
  }
  coefficients = t(as.data.frame(as.vector(coefficients)))
  colnames(coefficients)=col
  result = cbind(result, coefficients)
  #sigma
  sigma = summary$sigma
  result = cbind(result, sigma)
  #AICtab
  result = cbind(result, t(summary$AICtab))
  #append result into all_result
  if(cluster == 1){
    all_result = result
  } else { 
    #bind rows between all_result and result
    all_result = bind_rows(all_result,result)
  }
}

# row names of the results
row.names(all_result)= colnames(data)

# The first file is saved with column names
write.table(all_result,file=args[3],row.names = T,col.names = T)
