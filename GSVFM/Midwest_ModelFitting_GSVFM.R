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
source("GSVFM/app_funcs/agriculture_data_cv.R")

#### Loading the dataset ######
load("Data_Extraction_cleaning/RdataFiles/RawMidwestFundatPreProc_incl_irrig.RData")
load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcRegdat.RData")
load("Data_Extraction_cleaning/RdataFiles/MidwestScaledTriangulation.RData")

midw_years = length(unique(midw_regdat_new$Year))
NUM_ITERS = 9
fld_lim = 5

########################### Calculating the MSPE using percentages since now 1 score only explains 12% of the variation #######################
# coarse grid
var_perc = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
# fine-tuning
#var_perc = c(0.78,0.81,0.84,0.87,0.90,0.92,0.96,0.97,0.98)
deg = 3
midw_mspe_perc_var_trainMSE = list()
midw_mspe_perc_var_testMSPE = list()
num_scores_used = list()
for(i in 1:length(var_perc)){
  perc_var = var_perc[i]
  cat("The variance percentage is", perc_var, "\n")
  # Performing 3 repetitions for each
  midw_mspe_smooth_perc_var_trainMSE = data.frame(matrix(0, nr = 15, nc = 1))
  midw_mspe_smooth_perc_var_testMSPE = data.frame(matrix(0, nr = 15, nc = 5))
  colnames(midw_mspe_smooth_perc_var_trainMSE) = c("train_MSE")
  colnames(midw_mspe_smooth_perc_var_testMSPE) = c("test_MSPE", "test_MSPE_perc", "test_MSPE_var", "test_MSPE_mean_abs", "test_MSPE_mean_sqr")
  for (j in 1:3){
    ag.data.cv_otpt = ag.data_cv.fxn_nopresmooth(fundat_all = midw_fd_new, nonfd = midw_regdat_new,  n_yr = midw_years,
                                                 sp_tri = midw_scl_tri, n_fold = 5, iter = j, deg = deg, DE_MEAN_RESP = T, pred_vars = c("avgPRCP","irrig_prop"),
                                                 reval = c(0,1), fld_lim = fld_lim, thresh = perc_var)
    midw_mspe_smooth_perc_var_trainMSE[((j-1)*5 + 1):(5*j),] = ag.data.cv_otpt$train_MSE
    midw_mspe_smooth_perc_var_testMSPE[((j-1)*5 + 1):(5*j),] = ag.data.cv_otpt$test_MSPE
  }
  midw_mspe_perc_var_trainMSE[[i]] = midw_mspe_smooth_perc_var_trainMSE
  midw_mspe_perc_var_testMSPE[[i]] = midw_mspe_smooth_perc_var_testMSPE
  # The same number of scores is used for each repition so I only need to keep
  # track of it at the end
  num_scores_used[[i]] = ag.data.cv_otpt$num_used_harm
}
date = format(Sys.time(), "%m-%d-%Y")
file_name_trn = paste("OUT_Midwest_perc15to90by15_and95_TrainMSE_incl_irrig", date, ".RData", sep = "")
file_name_tst = paste("OUT_Midwest_perc15to90by15_and95_TestMSPE_incl_irrig", date, ".RData", sep = "")
file_name_numscores = paste("OUT_Midwest_perc15to90by15_and95_NumScores_incl_irrig", date, ".RData", sep = "")
# save(midw_mspe_perc_var_trainMSE, file = file_name_trn)
# save(midw_mspe_perc_var_testMSPE, file = file_name_tst)
# save(num_scores_used, file = file_name_numscores)
# Calcualting the averages across all folds
train_mse_avg = c()
test_mspe_avg = c()
for(i in 1:length(var_perc)){
  train_mse_avg[i] = mean(midw_mspe_perc_var_trainMSE[[i]][,1])
  test_mspe_avg[i] = mean(midw_mspe_perc_var_testMSPE[[i]]$test_MSPE)
}

# Saving the results across all percentage of variation values 
midw_test_mspe_cv = matrix(0, nrow = 15, nc = 2)
midw_test_mspe_cv[,1] = c(0.15,0.30,0.45,0.60,0.75,0.78,0.81,0.84,0.87,0.90,0.92,0.95, 0.96,0.97,0.98)
midw_test_mspe_cv[,2] = c(342.0221, 325.2699, 297.9933, 269.5976, 240.0083, 232.0209, 223.6730, 192.3272, 164.4852,
                          145.4474, 143.8232, 143.5626, 144.4571, 146.5286, 149.8723)
#save(midw_test_mspe_cv, file = "OUT_MidwestDataTestMSPE_CV_AllPercentages.RData")

############ Running 9 iterations using the optimal variance percentage ############
var_perc = c(0.90)
deg = 3
midw_mspe_perc_var_trainMSE = list()
midw_mspe_perc_var_testMSPE = list()
num_scores_used = list()
for(i in 1:length(var_perc)){
  perc_var = var_perc[i]
  cat("The variance percentage is", perc_var, "\n")
  # Performing 3 repetitions for each
  midw_mspe_smooth_perc_var_trainMSE = data.frame(matrix(0, nr = 15, nc = 1))
  midw_mspe_smooth_perc_var_testMSPE = data.frame(matrix(0, nr = 15, nc = 5))
  colnames(midw_mspe_smooth_perc_var_trainMSE) = c("train_MSE")
  colnames(midw_mspe_smooth_perc_var_testMSPE) = c("test_MSPE", "test_MSPE_perc", "test_MSPE_var", "test_MSPE_mean_abs", "test_MSPE_mean_sqr")
  for (j in 1:NUM_ITERS){
    ag.data.cv_otpt = ag.data_cv.fxn_nopresmooth(fundat_all = midw_fd_new, nonfd = midw_regdat_new,  n_yr = midw_years,
                                                 sp_tri = midw_scl_tri, n_fold = 5, iter = j, deg = deg, DE_MEAN_RESP = T, pred_vars = c("avgPRCP","irrig_prop"),
                                                 reval = c(0,1), fld_lim = fld_lim, thresh = perc_var)
    midw_mspe_smooth_perc_var_trainMSE[((j-1)*5 + 1):(5*j),] = ag.data.cv_otpt$train_MSE
    midw_mspe_smooth_perc_var_testMSPE[((j-1)*5 + 1):(5*j),] = ag.data.cv_otpt$test_MSPE
  }
  midw_mspe_perc_var_trainMSE[[i]] = midw_mspe_smooth_perc_var_trainMSE
  midw_mspe_perc_var_testMSPE[[i]] = midw_mspe_smooth_perc_var_testMSPE
  # The same number of scores is used for each repition so I only need to keep
  # track of it at the end
  num_scores_used[[i]] = ag.data.cv_otpt$num_used_harm
}
date = format(Sys.time(), "%m-%d-%Y")
file_name_trn = paste("OUT_Midwest_90perc_Opt_TrainMSE_incl_irrig", date, ".RData", sep = "")
file_name_tst = paste("OUT_Midwest_90perc_Opt_TestMSPE_incl_irrig", date, ".RData", sep = "")
file_name_numscores = paste("OUT_Midwest_90perc_Opt_NumScores_incl_irrig", date, ".RData", sep = "")
# save(midw_mspe_perc_var_trainMSE, file = file_name_trn)
# save(midw_mspe_perc_var_testMSPE, file = file_name_tst)
# save(num_scores_used, file = file_name_numscores)


#### Creating the cross-validation for var perc plot using ggplot #######
midw_train_mspe_cv = matrix(0, nrow = 15, nc = 2)
midw_test_mspe_cv = matrix(0, nrow = 15, nc = 2)
midw_train_mspe_cv[,1] = c(0.15,0.30,0.45,0.60,0.75,0.78,0.81,0.84,0.87,0.90,0.92,0.95, 0.96,0.97,0.98)
midw_train_mspe_cv[,2] = c(299.62, 288.12, 248.74, 204.32, 178.84, 155.17, 146.24, 117.69, 98.54, 85.10, 82.24, 77.49, 75.99,
                           74.24, 72.44)
midw_test_mspe_cv[,1] = c(0.15,0.30,0.45,0.60,0.75,0.78,0.81,0.84,0.87,0.90,0.92,0.95, 0.96,0.97,0.98)
midw_test_mspe_cv[,2] = c(342.0221, 325.2699, 297.9933, 269.5976, 240.0083, 232.0209, 223.6730, 192.3272, 164.4852,
                          145.4474, 143.8232, 143.5626, 144.4571, 146.5286, 149.8723)

# Organizing the data into a dataframe 
midw_mse_cv = data.frame(matrix(0, nrow = 30, ncol = 3))
midw_mse_cv[,1] = rep(c(0.15,0.30,0.45,0.60,0.75,0.78,0.81,0.84,0.87,0.90,0.92,0.95, 0.96,0.97,0.98))
midw_mse_cv[,2] = c(midw_train_mspe_cv[,2], midw_test_mspe_cv[,2])
midw_mse_cv[1:15,3] = "train"
midw_mse_cv[16:30,3] = "test"
colnames(midw_mse_cv) = c("perc_var", "MSE", "Set")
# Creating the plot
ggplot(data = midw_mse_cv, mapping = aes(x = perc_var, y = MSE, color = Set,linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("MSE")


















