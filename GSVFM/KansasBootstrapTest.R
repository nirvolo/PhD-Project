#source functions
source("common_utils/preprocessing_function.R") 
source('GSVFM/GSVCM_Kim_Wang/families.R')
source("GSVFM/GSVCM_Kim_Wang/gsvcm_est.R")
source("GSVFM/GSVCM_Kim_Wang/fit.gsvcm.R")
source("GSVFM/GSVCM_Kim_Wang/cv.gsvcm.R")
source("GSVFM/GSVCM_Kim_Wang/predict.gsvcm.R")
source("GSVFM/GSVCM_Kim_Wang/test.gsvcm.R")
source("GSVFM/util_funcs/train_fxn.R")
library(spatstat.geom) # this is for the function cross_dist.
library(MASS) # for mvrnorm
library("fda")
library(combinat)
library('BPST') # This package is for the basis function
library('MGLM')
library(plot3D)
library(devtools)
library(tidyverse)
library(ggplot2)
# Loading the Kansas data 
load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcFundat.RData")
load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcRegdat.RData")
load("Data_Extraction_cleaning/RdataFiles/KansasTriangulation.RData")

##### Setting up the predictor variables for the bootstrap hypothesis test #####
# Smoothing the entire kansas dataset. This will return me all the FPC scores for Kansas.
smoothed_kns_fd = fd_smooth_2steps_trn_tst(kns_fd_new, n_obs = nrow(kns_fd_new)/365, reval =  c(0,1), llam.lims = c(-15.5,5))
kns_comb_fd_lst = smoothed_kns_fd$comb_fd
### Centering the yield data by the mean by year ###
# Add column yld_mean_per_yr to train_nonfd
kns_nonfd_new$yld_mean_per_yr = ave(kns_nonfd_new$Yield, kns_nonfd_new$Year)
# add De-Meaned column to train
kns_nonfd_new$de_meaned_Yield = kns_nonfd_new$Yield - kns_nonfd_new$yld_mean_per_yr
# Fitting the training model
kns_mod = train_fxn_allfd(train_fd = kns_comb_fd_lst, train_nonfd = na.omit(kns_nonfd_new), full_fd = kns_comb_fd_lst,
                          thresh = 0.93, sp_tri = kns_tri2, use_full_fd4mean = F, d = 3, DE_MEAN_RESP = T, pred_vars = c("avgPRCP","irrig_prop"),
                          nharm = NULL)
kns_fpc_scores = kns_mod$fpc_scores
kns_locs = as.matrix(kns_nonfd_new[,c("long","lat")]) # the kns locations
B0_kns = basis(kns_tri2$V, kns_tri2$Tr, d = 3, r = 1, kns_locs)
Q2_kns = B0_kns$Q2
B_kns = B0_kns$B
# The Kansas data predictor variables which include the FPC scores
kns_pred_vars = data.matrix(cbind(rep(1,nrow(kns_nonfd_new)),kns_nonfd_new[,c("avgPRCP","irrig_prop")], kns_fpc_scores))
# Changing colnames for pred_vars
colnames(kns_pred_vars) = c("intercept", "avgPRCP", "irrig_prop", paste("PC",1), paste("PC",2), paste("PC",3), paste("PC",4),
                            paste("PC",5), paste("PC",6), paste("PC",7), paste("PC",8), paste("PC",9), paste("PC",10), paste("PC",111),
                            paste("PC",12), paste("PC",13), paste("PC",14), paste("PC",15), paste("PC",16))
y_resp = as.vector(kns_nonfd_new$de_meaned_Yield)

##### **Example running the bootstrap hypothesis test for Kansas data for one iteration** #######
d=3
r=1
lambda_start=0.00001
lambda_end=100
nlambda=10
lambda=exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))
test.result = test.gsvcm(y_resp, kns_pred_vars, kns_locs, kns_tri2$V, kns_tri2$Tr, d, r, lambda,
                          test_iter = 1, family = gaussian(), nB = 1)

################### Running Hypothesis Tests (Kansas Data) for all parameters with 500 reps ######################
d=3
r=1
lambda_start=0.00001
lambda_end=100
nlambda=10
lambda=exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))
n_reps = 500 # number of Bootstrap reps
# Creating a matrix to store the variable name and it's corresponding p-value
p.value.mat_kns = matrix(data = 0, nrow = ncol(kns_pred_vars), ncol = 3)
p.value.mat_kns[,1] = colnames(kns_pred_vars)
# To store the estimated statistic for each rep
test.stat_est = matrix(data = 0, nrow = 500, ncol = ncol(kns_pred_vars))
for(i in 1:nrow(p.value.mat_kns)){
  test.result = test.gsvcm(y_resp, kns_pred_vars, kns_locs, kns_tri2$V, kns_tri2$Tr, d, r, lambda,
                           test_iter = i, family = gaussian(), nB = n_reps)
  p.value.mat_kns[i,2] = test.result$obs.GQLR
  p.value.mat_kns[i,3] = test.result$pvalue
  test.stat_est[,i] = test.result$b.GQLR
}

# Saving the data files
date = format(Sys.time(), "%m-%d-%Y")
filename_pval = paste("P_Value_Matrix_Kansas", date, ".RData", sep = "")
filename_est.stat = paste("Est_HypTestStats_Matrix_Kansas", date, ".RData", sep = "")
#save(p.value.mat_kns, file = filename_pval)
#save(test.stat_est, file = filename_est.stat)









