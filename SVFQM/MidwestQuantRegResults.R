load("Data_Extraction_cleaning/Rdatafiles/RawMidwestFundatPreProc_incl_irrig.RData")
load("Data_Extraction_cleaning/Rdatafiles/RawMidwestRegdatPreProc_incl_irrig.RData")
load("Data_Extraction_cleaning/Rdatafiles/MidwestScaledTriangulation.RData")
source("common_utils/preprocessing_function.R")
source("SVFQM/app_funcs/agriculture_data_quantreg.R")

################################################################################################################
#################################### Midwest SVFQM CV for penalty parameter ######################################
################################################################################################################
lambda =  c(0.00025,0.0005,0.001,0.0015,0.003,0.005,0.01,0.025,0.05,0.075,0.10,0.15,0.20)
### For tau = 0.25
midw_svfqm_lambda_cv_tau0.25 = matrix(0, nrow = length(lambda), ncol = 3)
midw_svfqm_lambda_cv_tau0.25[,1] = lambda
colnames(midw_svfqm_lambda_cv_tau0.25) = c("lambda","train","test")
# this matrix is to store the QEP from all 5 folds
midw_svfqm_lam_QEP_tau0.25 = matrix(0, nrow = 11, ncol = length(lambda))
midw_svfqm_lam_QEP_tau0.25[1,] = lambda
rownames(midw_svfqm_lam_QEP_tau0.25) = c("lambda",rep("train",5),rep("test",5))
### For tau = 0.50
midw_svfqm_lambda_cv_tau0.5 = matrix(0, nrow = length(lambda), ncol = 3)
midw_svfqm_lambda_cv_tau0.5[,1] = lambda
colnames(midw_svfqm_lambda_cv_tau0.25) = c("lambda","train","test")
midw_svfqm_lam_QEP_tau0.5 = matrix(0, nrow = 11, ncol = length(lambda))
midw_svfqm_lam_QEP_tau0.5[1,] = lambda
rownames(midw_svfqm_lam_QEP_tau0.5) = c("lambda",rep("train",5),rep("test",5))
### For tau = 0.75
midw_svfqm_lambda_cv_tau0.75 = matrix(0, nrow = length(lambda), ncol = 3)
midw_svfqm_lambda_cv_tau0.75[,1] = lambda
colnames(midw_svfqm_lambda_cv_tau0.25) = c("lambda","train","test")
midw_svfqm_lam_QEP_tau0.75 = matrix(0, nrow = 11, ncol = length(lambda))
midw_svfqm_lam_QEP_tau0.75[1,] = lambda
rownames(midw_svfqm_lam_QEP_tau0.75) = c("lambda",rep("train",5),rep("test",5))

# using 60% FPCA variation for consistency across all quantiles
for(i in 1:length(lambda)){
  cat("The current value of lambda is",lambda[i],"\n")
  cat("   The current value of tau is", 0.25,"\n")
  midw_qr_res_pnlty_tau0.25 = ag.data_qr(fundat_all = midw_fd_new, nonfd = midw_regdat_new, thresh = 0.60, n_fold = 5, iter = 1, deg = 3, DE_MEAN_RESP = T,
                                        pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = 5, tau = 0.25, spat_triang = midw_scl_tri,
                                        comp_mod = F, svfqm = T, lambda_seq = c(lambda[i],1))
  cat("   The current value of tau is", 0.50,"\n")
  midw_qr_res_pnlty_tau0.5 = ag.data_qr(fundat_all = midw_fd_new, nonfd = midw_regdat_new, thresh = 0.60, n_fold = 5, iter = 1, deg = 3, DE_MEAN_RESP = T,
                                       pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = 5, tau = 0.5, spat_triang = midw_scl_tri,
                                       comp_mod = F, svfqm = T, lambda_seq = c(lambda[i],1))
  cat("   The current value of tau is", 0.75,"\n")
  midw_qr_res_pnlty_tau0.75 = ag.data_qr(fundat_all = midw_fd_new, nonfd = midw_regdat_new, thresh = 0.60, n_fold = 5, iter = 1, deg = 3, DE_MEAN_RESP = T,
                                        pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = 5, tau = 0.75, spat_triang = midw_scl_tri,
                                        comp_mod = F, svfqm = T, lambda_seq = c(lambda[i],1))
  ########## Storing the results #############
  midw_svfqm_lam_QEP_tau0.25[2:6,i] = midw_qr_res_pnlty_tau0.25$svqfm$train_QEP
  midw_svfqm_lam_QEP_tau0.25[7:11,i] = midw_qr_res_pnlty_tau0.25$svqfm$test_QEP
  ############
  midw_svfqm_lam_QEP_tau0.5[2:6,i] = midw_qr_res_pnlty_tau0.5$svqfm$train_QEP
  midw_svfqm_lam_QEP_tau0.5[7:11,i] = midw_qr_res_pnlty_tau0.5$svqfm$test_QEP
  ########
  midw_svfqm_lam_QEP_tau0.75[2:6,i] = midw_qr_res_pnlty_tau0.75$svqfm$train_QEP
  midw_svfqm_lam_QEP_tau0.75[7:11,i] = midw_qr_res_pnlty_tau0.75$svqfm$test_QEP
}

# Saving the results for each quantile
setwd("/Users/nirvoloshin/Documents/StatsResearch/Coding/QuantileRegression/MidwestResults")
date = format(Sys.time(), "%m-%d-%Y")
#Saving SVFQM results
file_name_tau0.25 = paste("MidwSVFQM_PnltyParamCV_tau0.25_", date, ".RData", sep = "")
file_name_tau0.50 = paste("MidwSVFQM_PnltyParamCV_tau0.50_", date, ".RData", sep = "")
file_name_tau0.75 = paste("MidwSVFQM_PnltyParamCV_tau0.75_", date, ".RData", sep = "")
# save(midw_svfqm_lam_QEP_tau0.25, file = file_name_tau0.25)
# save(midw_svfqm_lam_QEP_tau0.5, file = file_name_tau0.50)
# save(midw_svfqm_lam_QEP_tau0.75, file = file_name_tau0.75)

### Calculating the average QEP for each quantile and lambda value ####
# calculating the means
midw_svfqm_lambda_cv_tau0.25[,2] = apply(abs(midw_svfqm_lam_QEP_tau0.25[2:6,]),2,mean)
midw_svfqm_lambda_cv_tau0.25[,3] = apply(abs(midw_svfqm_lam_QEP_tau0.25[7:11,]),2,mean)
midw_svfqm_lambda_cv_tau0.5[,2] = apply(abs(midw_svfqm_lam_QEP_tau0.5[2:6,]),2,mean)
midw_svfqm_lambda_cv_tau0.5[,3] = apply(abs(midw_svfqm_lam_QEP_tau0.5[7:11,]),2,mean)
midw_svfqm_lambda_cv_tau0.75[,2] = apply(abs(midw_svfqm_lam_QEP_tau0.75[2:6,]),2,mean)
midw_svfqm_lambda_cv_tau0.75[,3] = apply(abs(midw_svfqm_lam_QEP_tau0.75[7:11,]),2,mean)

round(midw_svfqm_lambda_cv_tau0.25,digits = 5)
round(midw_svfqm_lambda_cv_tau0.5,digits = 5)
round(midw_svfqm_lambda_cv_tau0.75,digits = 5)


####### Creating a plot for the Lambda CV results #######
lambda = c(0.00025,0.0005,0.001,0.0015,0.003 ,0.005,0.01,0.025,0.05,0.075,0.10,0.15,0.20)
### 25th quantile
midw_tau0.25_svfqm_lambda_df = data.frame(matrix(0, nrow = 2*length(lambda), ncol = 3))
midw_tau0.25_svfqm_lambda_df[,1] = lambda
midw_tau0.25_svfqm_lambda_df[1:13,2] = midw_svfqm_lambda_cv_tau0.25[,2] # train values
midw_tau0.25_svfqm_lambda_df[1:13,3] = "train"
midw_tau0.25_svfqm_lambda_df[14:26,2] = midw_svfqm_lambda_cv_tau0.25[,3] # test values
midw_tau0.25_svfqm_lambda_df[14:26,3] = "test"
colnames(midw_tau0.25_svfqm_lambda_df) = c("lambda","QEP","Set")
### 50th quantile
midw_tau0.5_svfqm_lambda_df = data.frame(matrix(0, nrow = 2*length(lambda), ncol = 3))
midw_tau0.5_svfqm_lambda_df[,1] = lambda
midw_tau0.5_svfqm_lambda_df[1:13,2] = midw_svfqm_lambda_cv_tau0.5[,2] # train values
midw_tau0.5_svfqm_lambda_df[1:13,3] = "train"
midw_tau0.5_svfqm_lambda_df[14:26,2] = midw_svfqm_lambda_cv_tau0.5[,3] # test values
midw_tau0.5_svfqm_lambda_df[14:26,3] = "test"
colnames(midw_tau0.5_svfqm_lambda_df) = c("lambda","QEP","Set")
### 75th quantile
midw_tau0.75_svfqm_lambda_df = data.frame(matrix(0, nrow = 2*length(lambda), ncol = 3))
midw_tau0.75_svfqm_lambda_df[,1] = lambda
midw_tau0.75_svfqm_lambda_df[1:13,2] = midw_svfqm_lambda_cv_tau0.75[,2] # train values
midw_tau0.75_svfqm_lambda_df[1:13,3] = "train"
midw_tau0.75_svfqm_lambda_df[14:26,2] = midw_svfqm_lambda_cv_tau0.75[,3] # test values
midw_tau0.75_svfqm_lambda_df[14:26,3] = "test"
colnames(midw_tau0.75_svfqm_lambda_df) = c("lambda","QEP","Set")

ggplot(data = midw_tau0.25_svfqm_lambda_df, mapping = aes(x = lambda, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Nu") + ylab("QEP")

ggplot(data = midw_tau0.5_svfqm_lambda_df, mapping = aes(x = lambda, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Nu") + ylab("QEP")

ggplot(data = midw_tau0.75_svfqm_lambda_df, mapping = aes(x = lambda, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Nu") + ylab("QEP")

#############################################################################################################################
#################################### Midwest SVFQM and FQR CV for number of FPC scores ######################################
#############################################################################################################################
NUM_ITERS = 3 # only 3 iters since midw takes a long time to run
fld_lim = 5
ag_data_svfqm_QEP_res = list()
ag_data_fqr_QEP_res = list()
# coarse grid
perc_var = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
# fine-tuning changes depending on quantile
#perc_var = c(0.50,0.55,0.65,0.70)
tau = c(0.25,0.50,0.75)
# The lambda values for the SVFQM 
lambda_vals = c(0.025,0.01,0.0015)
tau_QEP_res_svfqm = list()
tau_QEP_res_fqr = list()
for(i in 1:length(tau)){
  lambda = lambda_vals[i]
  cat("The value of tau is ", tau[i], "\n")
  for(j in 1:length(perc_var)){
    cat("   The current percentage of variance is ", perc_var[j], "\n")
    # Organizing the matrices to store the train and test QEP for each percentage
    ag_data_QEP_svfqm = data.frame(matrix(0, nr = 15, nc = 2))
    colnames(ag_data_QEP_svfqm) = c("train QEP", "test QEP")
    ag_data_QEP_fqr = data.frame(matrix(0, nr = 15, nc = 2))
    colnames(ag_data_QEP_fqr) = c("train QEP", "test QEP")
    for (k in 1:NUM_ITERS){
      midw_qr_res_otpt = ag.data_qr(fundat_all = midw_fd_new, nonfd = midw_regdat_new, thresh = perc_var[j], n_fold = 5, iter = k, deg = 3, DE_MEAN_RESP = T,
                                   pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = fld_lim, tau = tau[i], spat_triang = midw_scl_tri,
                                   comp_mod = T, svfqm = T, lambda_seq = c(lambda,2))
      # Updating SVFQM results
      ag_data_QEP_svfqm[((k-1)*5 + 1):(5*k),1] = midw_qr_res_otpt$svqfm$train_QEP
      ag_data_QEP_svfqm[((k-1)*5 + 1):(5*k),2] = midw_qr_res_otpt$svqfm$test_QEP
      # Updating FQR results
      ag_data_QEP_fqr[((k-1)*5 + 1):(5*k),1] = midw_qr_res_otpt$fqr$train_MSE
      ag_data_QEP_fqr[((k-1)*5 + 1):(5*k),2] = midw_qr_res_otpt$fqr$test_MSE
    }
    ag_data_svfqm_QEP_res[[j]] = ag_data_QEP_svfqm
    ag_data_fqr_QEP_res[[j]] = ag_data_QEP_fqr
    # The same number of scores is used for each repition so I only need to keep
    # track of it at the end
    num_scores_used_svfqm[[j]] = midw_qr_res_otpt$svqfm$num_used_harm
    num_scores_used_fqr[[j]] = midw_qr_res_otpt$fqr$num_used_harm
  }
  # Store results for the tau'th quantile
  tau_QEP_res_svfqm[[i]] = ag_data_svfqm_QEP_res
  tau_QEP_res_fqr[[i]] = ag_data_fqr_QEP_res
}

# Saving the SVFQM and FQR results
date = format(Sys.time(), "%m-%d-%Y")
# Saving SVFQM results
file_name_svfqm = paste("UpdtLambda_Midw_SVFQM_CV_W_Pnlty_Res_tau0.25_0.50_0.75_PercVar0.15to0.90_0.95", date, ".RData", sep = "")
#save(tau_QEP_res_svfqm, file = file_name_svfqm)
file_name_svfqm_num_scores = paste("OUT_Midw_SVFQM_NumScores_0.80_85", date, ".RData", QEP = "")
#save(num_scores_used_svfqm, file = file_name_svfqm_num_scores)
file_name_fqr = paste("OUT_Midw_UpdatedFQR_CV_W_Pnlty_Res_tau0.25_0.50_0.75_PercVar0.15-0.90_0.95", date, ".RData", QEP = "")
#save(tau_QEP_res_fqr, file = file_name_fqr)
file_name_fqr_num_scores = paste("OUT_Midw_FQR_NumScores_0.15-0.90_0.95", date, ".RData", QEP = "")
#save(num_scores_used_fqr, file = file_name_fqr_num_scores)


####### Calculating the average QEP for the percentiles and quantiles above ########
perc_var = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
tau_0.25_avg_QEP_svfqm = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.25_avg_QEP_svfqm[,1] = perc_var
colnames(tau_0.25_avg_QEP_svfqm) = c("perc_var","train_QEP", "test_QEP")
tau_0.25_avg_QEP_fqr = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.25_avg_QEP_fqr[,1] = perc_var
colnames(tau_0.25_avg_QEP_fqr) = c("perc_var","train_QEP", "test_QEP")
tau_0.5_avg_QEP_svfqm = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.5_avg_QEP_svfqm[,1] = perc_var
colnames(tau_0.5_avg_QEP_svfqm) = c("perc_var","train_QEP", "test_QEP")
tau_0.5_avg_QEP_fqr = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.5_avg_QEP_fqr[,1] = perc_var
colnames(tau_0.5_avg_QEP_fqr) = c("perc_var","train_QEP", "test_QEP")
tau_0.75_avg_QEP_svfqm = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.75_avg_QEP_svfqm[,1] = perc_var
colnames(tau_0.5_avg_QEP_svfqm) = c("perc_var","train_QEP", "test_QEP")
tau_0.75_avg_QEP_fqr = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.75_avg_QEP_fqr[,1] = perc_var
colnames(tau_0.75_avg_QEP_fqr) = c("perc_var","train_QEP", "test_QEP")

for(j in 1:length(perc_var)){
  tau_0.25_avg_QEP_svfqm[j,2] = mean(abs(tau_QEP_res_svfqm_tau0.25[[1]][[j]]$`train QEP`))
  tau_0.25_avg_QEP_svfqm[j,3] = mean(abs(tau_QEP_res_svfqm_tau0.25[[1]][[j]]$`test QEP`))
  tau_0.25_avg_QEP_fqr[j,2] = mean(abs(tau_QEP_res_fqr[[1]][[j]]$`train QEP`))
  tau_0.25_avg_QEP_fqr[j,3] = mean(abs(tau_QEP_res_fqr[[1]][[j]]$`test QEP`))                                
  
  tau_0.5_avg_QEP_svfqm[j,2] = mean(abs(tau_QEP_res_svfqm_tau0.50[[2]][[j]]$`train QEP`))
  tau_0.5_avg_QEP_svfqm[j,3] = mean(abs(tau_QEP_res_svfqm_tau0.50[[2]][[j]]$`test QEP`))
  tau_0.5_avg_QEP_fqr[j,2] = mean(abs(tau_QEP_res_fqr[[2]][[j]]$`train QEP`))
  tau_0.5_avg_QEP_fqr[j,3] = mean(abs(tau_QEP_res_fqr[[2]][[j]]$`test QEP`))    
  
  tau_0.75_avg_QEP_svfqm[j,2] = mean(abs(tau_QEP_res_svfqm_tau0.75[[3]][[j]]$`train QEP`))
  tau_0.75_avg_QEP_svfqm[j,3] = mean(abs(tau_QEP_res_svfqm_tau0.75[[3]][[j]]$`test QEP`))
  tau_0.75_avg_QEP_fqr[j,2] = mean(abs(tau_QEP_res_fqr[[3]][[j]]$`train QEP`))
  tau_0.75_avg_QEP_fqr[j,3] = mean(abs(tau_QEP_res_fqr[[3]][[j]]$`test QEP`)) 
  
}

round(tau_0.25_avg_QEP_svfqm,digits = 5)
round(tau_0.5_avg_QEP_svfqm,digits = 5)
round(tau_0.75_avg_QEP_svfqm, digits = 5)

round(tau_0.25_avg_QEP_fqr,digits = 5)
round(tau_0.5_avg_QEP_fqr,digits = 5)
round(tau_0.75_avg_QEP_fqr, digits = 5)

#################################### Midwest SVFQM and FQR CV Perc Var Plots ####################################
###### Results for tau = 0.25 ######
#save(tau_0.25_avg_QEP_svfqm_new, file = "MidwSVFQM_tau0.25_NumScoresCV_Matrix.RData")
load("Data_Extraction_cleaning/RDataFiles/SVFQM_Results/MidwSVFQM_tau0.25_NumScoresCV_Matrix.RData")
midw_QEP_scores_cv_svfqm_tau0.25 = matrix(0, nrow = 22, 2)
midw_QEP_scores_cv_svfqm_tau0.25[1:11,1] = tau_0.25_avg_QEP_svfqm_new[,1]
midw_QEP_scores_cv_svfqm_tau0.25[12:22,1] = tau_0.25_avg_QEP_svfqm_new[,1]
midw_QEP_scores_cv_svfqm_tau0.25[1:11,2] = tau_0.25_avg_QEP_svfqm_new[,2] # the train QEP
midw_QEP_scores_cv_svfqm_tau0.25[12:22,2] = tau_0.25_avg_QEP_svfqm_new[,3]  # test QEP
midw_QEP_scores_cv_svfqm_tau0.25_df = data.frame(midw_QEP_scores_cv_svfqm_tau0.25)
midw_QEP_scores_cv_svfqm_tau0.25_df[1:11,"Set"] = "train"
midw_QEP_scores_cv_svfqm_tau0.25_df[12:22,"Set"] = "test"
colnames(midw_QEP_scores_cv_svfqm_tau0.25_df) = c("PercVar","QEP", "Set")
ggplot(data = midw_QEP_scores_cv_svfqm_tau0.25_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")

###### Results for tau = 0.50 ######
##### SVFQM results
#save(tau_0.5_avg_QEP_svfqm, file = "MidwSVFQM_tau0.50_NumScoresCV_Matrix.RData")
load("Data_Extraction_cleaning/RDataFiles/SVFQM_Results/MidwSVFQM_tau0.25_NumScoresCV_Matrix.RData")
midw_QEP_scores_cv_svfqm_tau0.50 = matrix(0, nrow = 14, 2)
midw_QEP_scores_cv_svfqm_tau0.50[1:7,1] = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
midw_QEP_scores_cv_svfqm_tau0.50[8:14,1] = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
midw_QEP_scores_cv_svfqm_tau0.50[1:7,2] = tau_0.5_avg_QEP_svfqm[,2] # the train QEP
midw_QEP_scores_cv_svfqm_tau0.50[8:14,2] = tau_0.5_avg_QEP_svfqm[,3]  # test QEP
midw_QEP_scores_cv_svfqm_tau0.50_df = data.frame(midw_QEP_scores_cv_svfqm_tau0.50)
midw_QEP_scores_cv_svfqm_tau0.50_df[1:7,"Set"] = "train"
midw_QEP_scores_cv_svfqm_tau0.50_df[8:14,"Set"] = "test"
colnames(midw_QEP_scores_cv_svfqm_tau0.50_df) = c("PercVar","QEP", "Set")
ggplot(data = midw_QEP_scores_cv_svfqm_tau0.50_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")

###### Results for tau = 0.75 ######
#### SVFQM results
tau_0.75_avg_QEP_svfqm_new = rbind(tau_0.75_avg_QEP_svfqm[1:3,],tau_0.75_avg_QEP_svfqm_fine.tune[1:2,],
                                   tau_0.75_avg_QEP_svfqm[4,],tau_0.75_avg_QEP_svfqm_fine.tune[3:4,],
                                   tau_0.75_avg_QEP_svfqm[5:7,])
#save(tau_0.75_avg_QEP_svfqm_new,file = "MidwSVFQM_tau0.75_NumScoresCV_Matrix.RData")
load("Data_Extraction_cleaning/RDataFiles/SVFQM_Results/MidwSVFQM_tau0.25_NumScoresCV_Matrix.RData")
midw_QEP_scores_cv_svfqm_tau0.75 = matrix(0, nrow = 22, 2)
midw_QEP_scores_cv_svfqm_tau0.75[1:11,1] = tau_0.75_avg_QEP_svfqm_new[,1]
midw_QEP_scores_cv_svfqm_tau0.75[12:22,1] = tau_0.75_avg_QEP_svfqm_new[,1]
midw_QEP_scores_cv_svfqm_tau0.75[1:11,2] = tau_0.75_avg_QEP_svfqm_new[,2] # the train QEP
midw_QEP_scores_cv_svfqm_tau0.75[12:22,2] = tau_0.75_avg_QEP_svfqm_new[,3]  # test QEP
midw_QEP_scores_cv_svfqm_tau0.75_df = data.frame(midw_QEP_scores_cv_svfqm_tau0.75)
midw_QEP_scores_cv_svfqm_tau0.75_df[1:11,"Set"] = "train"
midw_QEP_scores_cv_svfqm_tau0.75_df[12:22,"Set"] = "test"
colnames(midw_QEP_scores_cv_svfqm_tau0.75_df) = c("PercVar","QEP", "Set")
ggplot(data = midw_QEP_scores_cv_svfqm_tau0.75_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")


#################################### Final Results for Midwest Data SVFQM 9 ITERS ######################################
NUM_ITERS = 9
fld_lim = 5
# Creating a matrix with the quantile and percent variation
# for each model and quantile. For example, 1st row, corresponds
# to tau = 0.25, SVFQM and optimal percentage
midw_smooth_params = matrix(0, nrow = 6, ncol = 3)
midw_smooth_params[,1] = c(0.25,0.25,0.50,0.50,0.75,0.75) # the quantile value
midw_smooth_params[,2] = c(0.60,0.15,0.95,0.80,0.55,0.90) # optimal variance percentage
midw_smooth_params[,3] = c(0.025,0,0.01,0,0.0015,0) # optimal lambda
midw_mods = c("SVFQM","FQR","SVFQM","FQR","SVFQM","FQR")
midw_final_res = list()
for(i in 1:nrow(midw_smooth_params)){
  tau = midw_smooth_params[i,1]
  perc_var = midw_smooth_params[i,2]
  mod = midw_mods[i]
  lambda = midw_smooth_params[i,3]
  cat("The current model is", mod, "\n")
  cat("  The value of tau is", tau, "\n")
  cat("    The current percentage of variance is", perc_var, "\n")
  # Organizing the matrices to store the train and test MSE for each percentage
  ag_data_mse = data.frame(matrix(0, nr = 45, nc = 2))
  colnames(ag_data_mse) = c("train QEP", "test QEP")
  if(mod == "SVFQM"){
    comp_mod = F
    svfqm = T
  } else{
    comp_mod = T
    svfqm = F
  }
  
  for(k in 1:NUM_ITERS){
    midw_qr_res_otpt = ag.data_qr(fundat_all = midw_fd_new, nonfd = midw_regdat_new, thresh = perc_var, n_fold = 5, iter = k, deg = 3,
                                  DE_MEAN_RESP = T, pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = fld_lim, tau = tau,
                                  spat_triang = midw_scl_tri, comp_mod = comp_mod, svfqm = svfqm, lambda_seq = c(lambda,5))
  
    if(svfqm){
      # Updating SVFQM results
      ag_data_mse[((k-1)*5 + 1):(5*k),1] = midw_qr_res_otpt$svqfm$train_QEP
      ag_data_mse[((k-1)*5 + 1):(5*k),2] = midw_qr_res_otpt$svqfm$test_QEP
    } else{
      # Updating FQR results
      ag_data_mse[((k-1)*5 + 1):(5*k),1] = midw_qr_res_otpt$fqr$train_QEP
      ag_data_mse[((k-1)*5 + 1):(5*k),2] = midw_qr_res_otpt$fqr$test_QEP
    }
  }
  # Updating the results for the corresponding combination of model, tau, and perc var
  midw_final_res[[i]] = ag_data_mse
}

date = format(Sys.time(), "%m-%d-%Y")
file_name = paste("OUT_Midw_SVFQM_FQR_W_Pnlty_tau0.25_0.40_0.45_0.75_OptScores", date, ".RData", QEP = "")
#save(midw_final_res, file = file_name)

#### Calculating the final results ####
midw_smooth_params_res = matrix(0, nrow = 6, ncol = 5)
midw_smooth_params_res[,1] = c(0.25,0.25,0.50,0.50,0.75,0.75) # the quantile value
midw_smooth_params_res[,2] = c(0.60,0.15,0.95,0.80,0.55,0.90) # optimal variance percentage
midw_smooth_params_res[,3] = c("SVFQM", "FQR", "SVFQM", "FQR", "SVFQM", "FQR")
colnames(midw_smooth_params_res) = c("quantile","perc_var","model","train QEP", "test QEP")

# Calculating the averages
for(i in 1:length(midw_final_res)){
  # Calculating avg train and test QEP
  midw_smooth_params_res[i,4] = mean(abs(midw_final_res[[i]]$`train QEP`))
  midw_smooth_params_res[i,5] = mean(abs(midw_final_res[[i]]$`test QEP`))
}







