source("GSVFM/sim_funcs/GSVFM_SimSettings.R")
source("GSVFM/app_funcs/agriculture_data_cv.R")

##################### Subsetting the Kansas data from the fd and nonfd midwest data for years 1999-2011 ##################### 
load("Data_extraction_cleaning/Rdatafiles/RawMidwestFundatBeforePreProcFxn_incl_irrig.RData") # midw_fd
load("Data_extraction_cleaning/Rdatafiles/RawMidwestRegdatBeforePreProcFxn_incl_irrig.RData") # midw_nonfd
kns_fd = midw_fd[midw_fd$State == "kansas" & midw_fd$Year %in% 1999:2011,]
kns_nonfd = midw_regdat[midw_regdat$State == "kansas" & midw_regdat$Year %in% 1999:2011,]
preproces_otpt_new_kns = preprocessing_nosmooth(fun_dat = kns_fd, non_fd = kns_nonfd, states = c("Kansas"), vars = c("Yield", "Area", 
                                                        "Irrigated.actualprop", "avgPRCP", "CountyI", "Year", "State"), st_bnd = list(kns_mat))
kns_fd_new = preproces_otpt_new_kns$fd_dat
kns_nonfd_new = preproces_otpt_new_kns$non_fd
# Adding observation Id to fd and nonfd
kns_nonfd_new[,"Id"] = seq(1, nrow(kns_nonfd_new))
kns_fd_new[,"Id"] = NA
for(i in 1:(nrow(kns_fd_new)/365)){
  kns_fd_new[((i-1)*365 + 1):(365*i),"Id"] = i
}
colnames(kns_fd_new) = c("State", "County", "Year",  "DateI",  "TMAX", "TMIN", "CountyI", "county", "long", "lat",
                          "Id")
# Saving the preprocessed Kansas data
# save(kns_nonfd_new, file = "NewKansasFromMidw_PreProcRegdat.RData")
# save(kns_fd_new, file = "NewKansasFromMidw_PreProcFundat.RData")

####### **Example of running real data function** #######
kns_data_gsvfm_ex = ag.data_cv.fxn_nopresmooth(fundat_all = kns_fd_new, nonfd = kns_nonfd_new, thresh = 0.50, n_fold = 5, iter = 1, 
                                               deg = 3, DE_MEAN_RESP = T, pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = 5)

######################################### Cross-Validation for Number of FPC Scores with Kansas data #####################################
NUM_ITERS = 9
fld_lim = 5
ag_data_mspe_list = list()
ag_data_train_mse_list = list()
# coarse search of optimal FPC scores to include
# perc_var = c(0.15,0.3,0.45,0.6,0.75,0.90,0.95)
# Fine tune after results from above coarse search
perc_var = c(0.93,0.94)
for(i in 1:length(perc_var)){
  cat("The current percentage of variance is ", perc_var[i], "\n") 
  ag_data_mspe_percvar = data.frame(matrix(0, nr = 45, nc = 5))
  ag_data_mse_perc_var = data.frame(matrix(0, nr = 45, nc = 5))
  colnames(ag_data_mspe_percvar) = c("test_MSE", "test_mse_perc", "test_mse_var", "test_mse_mean_abs", "test_mse_mean_sqr")
  colnames(ag_data_mse_perc_var) = c("train_MSE")
  for (j in 1:NUM_ITERS){
    ag.data.cv_otpt = ag.data_cv.fxn_nopresmooth(fundat_all = kns_fd_new, nonfd = kns_nonfd_new, thresh = perc_var[i], n_fold = 5, iter = j, 
                                                 deg = 3, DE_MEAN_RESP = T, pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = fld_lim)
    ag_data_mse_perc_var[((j-1)*5 + 1):(5*j),] = ag.data.cv_otpt$train_MSE
    ag_data_mspe_percvar[((j-1)*5 + 1):(5*j),] = ag.data.cv_otpt$test_MSPE
  }
  ag_data_train_mse_list[[i]] = ag_data_mse_perc_var
  ag_data_mspe_list[[i]] = ag_data_mspe_percvar
}

## Saving the CV results for different percentages ##
# save(ag_data_train_mse_list,file = "NewKansasData_CV_15to95percent_TrainMSE.RData")
# save(ag_data_mspe_list,file = "NewKansasData_CV_15to95percent_TestMSPE.RData")

# save(ag_data_train_mse_list,file = "NewKansasData_CV_96to98percent_TrainMSE.RData")
# save(ag_data_mspe_list,file = "NewKansasData_CV_96to98percent_TestMSPE.RData")

# save(ag_data_train_mse_list,file = "NewKansasData_CV_78to92percent_TrainMSE.RData")
# save(ag_data_mspe_list,file = "NewKansasData_CV_78to92percent_TestMSPE.RData")

# save(ag_data_train_mse_list,file = "NewKansasData_CV_93to94percent_TrainMSE.RData")
# save(ag_data_mspe_list,file = "NewKansasData_CV_93to94percent_TestMSPE.RData")


## Calculating the averages for each percentage ##
kns_train_mse = c()
kns_test_mspe = c()
for(i in 1:length(perc_var)){
  kns_train_mse[i] = mean(ag_data_train_mse_list[[i]][,1])
  kns_test_mspe[i] = mean(ag_data_mspe_list[[i]]$test_MSE)
}

####### Summarizing all the percentages in a plot #############
# The training MSE for all percentages
kns_train_mse_scores_cv = matrix(0, nrow = 17, 2)
kns_train_mse_scores_cv[,1] = c(0.15,0.30,0.45,0.60,0.75, 0.78,0.81,0.84,0.87,0.90,0.92,0.93,
                                0.94, 0.95, 0.96,0.97,0.98)
kns_train_mse_scores_cv[,2] = c(238.8364, 233.3671, 239.6824, 235.6421, 188.9131, 188.9131, 178.7256, 
                                168.5767, 156.2550, 155.4023, 105.1596, 103.0284, 104.1101, 103.0866, 109.3028, 107.6052, 102.1385)
# The test MSE for all percentages
kns_test_mspe_scores_cv = matrix(0, nrow = 17, 2)
kns_test_mspe_scores_cv[,1] = c(0.15,0.30,0.45,0.60,0.75,0.78,0.81,0.84,0.87,0.90,0.92,0.93,0.94, 0.95, 0.96,0.97,0.98)
kns_test_mspe_scores_cv[,2] = c(322.0854,331.1668,342.0714, 323.3530, 286.1240, 286.1240, 276.2257, 271.1997, 258.3294, 
                                257.3755, 196.0433, 195.5272, 196.7224, 198.6950, 199.7254, 206.3242, 214.3263)
plot(kns_test_mspe_scores_cv[,1], kns_test_mspe_scores_cv[,2], type = "b", pch = 20, cex = 0.8, main = "Test MSPE vs. Variance Percentage",
     xlab = "Variance Percentage", ylab = "MSPE")

# Creating a plot summarizing all the MSE results for Kansas
library(tidyverse)
# Matrix to store all the data
kns_mse_scores_cv = data.frame(matrix(c(kns_train_mse_scores_cv[,1], kns_test_mspe_scores_cv[,1], kns_train_mse_scores_cv[,2], kns_test_mspe_scores_cv[,2]), nrow = 34, ncol = 2))
kns_mse_scores_cv[1:17,"Set"] = "train"
kns_mse_scores_cv[18:34,"Set"] = "test"
colnames(kns_mse_scores_cv) = c("PercVar","MSE", "Set")
ggplot(data = kns_mse_scores_cv, mapping = aes(x = PercVar, y = MSE, color = Set,linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("MSE")

# Saving the file with all the CV results
#save(kns_test_mspe_scores_cv, file = "OUT_KansasDataTestMSPE_CV_AllPercentages.RData")




