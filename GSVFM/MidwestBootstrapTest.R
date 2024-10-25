#source functions
source("common_utils/preprocessing_function.R") 
source('GSVFM/GSVCM_Kim_Wang/families.R')
source("GSVFM/GSVCM_Kim_Wang/gsvcm_est.R")
source("GSVFM/GSVCM_Kim_Wang/fit.gsvcm.R")
source("GSVFM/GSVCM_Kim_Wang/cv.gsvcm.R")
source("GSVFM/GSVCM_Kim_Wang/predict.gsvcm.R")
source("GSVFM/GSVCM_Kim_Wang/test.gsvcm.R")
source("GSVFM/util_funcs/train_fxn.R")
load("Data_Extraction_cleaning/RdataFiles/RawMidwestFundatPreProc_incl_irrig.RData") # name in R is "midw_fd_new"
load("Data_Extraction_cleaning/RdataFiles/RawMidwestRegdatPreProc_incl_irrig.RData") # name in R is "midw_regdat_new"
load("Data_Extraction_cleaning/RdataFiles/MidwestScaledTriangulation.RData")

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

##### Setting up the predictor variables for the bootstrap hypothesis test #####
# Smoothing the entire midwest dataset. This will return me all the FPC scores for Kansas.
smoothed_midw_fd = fd_smooth_2steps_trn_tst(midw_fd_new, n_obs = nrow(midw_fd_new)/365, reval =  c(0,1), llam.lims = c(-15.5,5))
midw_comb_fd_lst = smoothed_midw_fd$comb_fd

### Centering the yield data by the mean by year ###
# Add column yld_mean_per_yr to train_nonfd
midw_regdat_new$yld_mean_per_yr = ave(midw_regdat_new$Yield, midw_regdat_new$Year)
# add De-Meaned column to train
midw_regdat_new$de_meaned_Yield = midw_regdat_new$Yield - midw_regdat_new$yld_mean_per_yr
# Fitting the training. Using 90% variation since that's what I found to be optimal for midwest data
midw_mod = train_fxn_allfd(train_fd = midw_comb_fd_lst, train_nonfd = na.omit(midw_regdat_new), full_fd = midw_comb_fd_lst,
                          thresh = 0.90, sp_tri = midw_scl_tri, use_full_fd4mean = F, d = 3, DE_MEAN_RESP = T, pred_vars = c("avgPRCP","irrig_prop"),
                          nharm = NULL)
midw_fpc_scores = midw_mod$fpc_scores

midw_locs = as.matrix(midw_regdat_new[,c("long","lat")]) # the midw locations
B0_midw = basis(midw_scl_tri$V, midw_scl_tri$Tr, d = 3, r = 1, midw_locs)
Q2_midw = B0_midw$Q2
B_midw = B0_midw$B
# The Kansas data predictor variables which include the FPC scores
midw_pred_vars = data.matrix(cbind(rep(1,nrow(midw_regdat_new)),midw_regdat_new[,c("avgPRCP","irrig_prop")], midw_fpc_scores))
# Changing colnames for pred_vars
colnames(midw_pred_vars) = c("intercept", "avgPRCP", "irrig_prop", paste("PC",1), paste("PC",2), paste("PC",3), paste("PC",4),
                            paste("PC",5), paste("PC",6), paste("PC",7), paste("PC",8), paste("PC",9), paste("PC",10), paste("PC",111),
                            paste("PC",12), paste("PC",13), paste("PC",14), paste("PC",15), paste("PC",16), paste("PC",17), paste("PC",18),
                            paste("PC",19), paste("PC",20), paste("PC",21))
y_resp = as.vector(midw_regdat_new$de_meaned_Yield)

##### **Example running the bootstrap hypothesis test for Kansas data for one iteration** #######
d=3
r=1
lambda_start=0.00001
lambda_end=100
nlambda=10
lambda=exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))
test.result = test.gsvcm(y_resp, midw_pred_vars, midw_locs,midw_scl_tri$V, midw_scl_tri$Tr, d, r, lambda,
                         test_iter = 1, family = gaussian(), nB = 1)

################### Running Hypothesis Tests (Midwest Data) for parameters 1-5 on the server with 500 reps ######################
d=3
r=1
lambda_start=0.00001
lambda_end=100
nlambda=10
lambda=exp(seq(log(lambda_start),log(lambda_end),length.out=nlambda))
n_reps = 500 # number of Bootstrap reps
## **The bootstrap tests can be split up into different files to speed up computation** ##
n_vars = ncol(midw_pred_vars)
# Creating a matrix to store the variable name and it's corresponding p-value
# Creating four rows since I only want to do the first 4 variables
p.value.mat_midw = matrix(data = 0, nrow = n_vars, ncol = 3)
p.value.mat_midw[,1] = colnames(midw_pred_vars)[1:5]
# # To store the estimated statistic for each rep. Only doing this for 5 variables
test.stat_est = matrix(data = 0, nrow = n_reps, ncol = n_vars)
for(i in 1:n_vars){
    test.result = test.gsvcm(y_resp, midw_pred_vars, midw_locs, midw_scl_tri$V, midw_scl_tri$Tr, d, r, lambda,
                         test_iter = i, family = gaussian(), nB = n_reps)
    p.value.mat_midw[i,2] = test.result$obs.GQLR
    p.value.mat_midw[i,3] = test.result$pvalue
    test.stat_est[,i] = test.result$b.GQLR
}

##########################################################################################
# The hypothesis tests were originally split into separate files where each file contained 
# 4-5 variables. The code is the same but the the for loop over n_vars changes. The file
# saved here is for variables 1 to 5 but the rest of the results are located in 
# Data_Extraction_cleaning/RdataFiles/BootstrapHypothesisTests
#########################################################################################
# # Saving the data files
date = format(Sys.time(), "%m-%d-%Y")
filename_pval = paste("P_Value_Matrix_Midw_Vars1to5", date, ".RData", sep = "")
filename_est.stat = paste("Est_HypTestStats_Matrix_Midw_Vars1to5", date, ".RData", sep = "")
# save(p.value.mat_midw, file = filename_pval)
# save(test.stat_est, file = filename_est.stat)


















