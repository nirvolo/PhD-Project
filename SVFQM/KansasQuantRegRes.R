load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcFundat.RData")
load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcRegdat.RData")
load("Data_Extraction_cleaning/RdataFiles/KansasTriangulation.RData")
source("SVFQM/app_funcs/agriculture_data_quantreg.R")

###### Example of running the ag.data_qr function with the Kansas data ########
kns_qr_res = ag.data_qr(fundat_all = kns_fd_new, nonfd = kns_nonfd_new, thresh = 0.60, n_fold = 5, iter = 1, deg = 3, DE_MEAN_RESP = T,
                        pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = 1, tau = 0.25, spat_triang = kns_tri2,
                        comp_mod = T, svfqm = T)

#################################################################################################################
#################################### Kansas SVFQM CV for penalty parameter ######################################
#################################################################################################################
lambda =  c(0.00025,0.0005,0.001,0.0015,0.003 ,0.005,0.01,0.025,0.05,0.075,0.10,0.15,0.20)
### For tau = 0.25
kns_svfqm_lambda_cv_tau0.25 = matrix(0, nrow = length(lambda), ncol = 3)
kns_svfqm_lambda_cv_tau0.25[,1] = lambda
colnames(kns_svfqm_lambda_cv_tau0.25) = c("lambda","train","test")
# this matrix is to store the QEP from all 5 folds
kns_svfqm_lam_QEP_tau0.25 = matrix(0, nrow = 11, ncol = length(lambda))
kns_svfqm_lam_QEP_tau0.25[1,] = lambda
rownames(kns_svfqm_lam_QEP_tau0.25) = c("lambda",rep("train",5),rep("test",5))
### For tau = 0.50
kns_svfqm_lambda_cv_tau0.5 = matrix(0, nrow = length(lambda), ncol = 3)
kns_svfqm_lambda_cv_tau0.5[,1] = lambda
colnames(kns_svfqm_lambda_cv_tau0.25) = c("lambda","train","test")
kns_svfqm_lam_QEP_tau0.5 = matrix(0, nrow = 11, ncol = length(lambda))
kns_svfqm_lam_QEP_tau0.5[1,] = lambda
rownames(kns_svfqm_lam_QEP_tau0.5) = c("lambda",rep("train",5),rep("test",5))
### For tau = 0.75
kns_svfqm_lambda_cv_tau0.75 = matrix(0, nrow = length(lambda), ncol = 3)
kns_svfqm_lambda_cv_tau0.75[,1] = lambda
colnames(kns_svfqm_lambda_cv_tau0.25) = c("lambda","train","test")
kns_svfqm_lam_QEP_tau0.75 = matrix(0, nrow = 11, ncol = length(lambda))
kns_svfqm_lam_QEP_tau0.75[1,] = lambda
rownames(kns_svfqm_lam_QEP_tau0.75) = c("lambda",rep("train",5),rep("test",5))


for(i in 1:length(lambda)){
  cat("The current value of lambda is",lambda[i],"\n")
  cat("   The current value of tau is", 0.25,"\n")
  kns_qr_res_pnlty_tau0.25 = ag.data_qr(fundat_all = kns_fd_new, nonfd = kns_nonfd_new, thresh = 0.60, n_fold = 5, iter = 1, deg = 3, DE_MEAN_RESP = T,
                                        pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = 5, tau = 0.25, spat_triang = kns_tri2,
                                        comp_mod = F, svfqm = T, lambda_seq = c(lambda[i],1))
  cat("   The current value of tau is", 0.50,"\n")
  kns_qr_res_pnlty_tau0.5 = ag.data_qr(fundat_all = kns_fd_new, nonfd = kns_nonfd_new, thresh = 0.60, n_fold = 5, iter = 1, deg = 3, DE_MEAN_RESP = T,
                                        pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = 5, tau = 0.5, spat_triang = kns_tri2,
                                        comp_mod = F, svfqm = T, lambda_seq = c(lambda[i],1))
  cat("   The current value of tau is", 0.75,"\n")
  kns_qr_res_pnlty_tau0.75 = ag.data_qr(fundat_all = kns_fd_new, nonfd = kns_nonfd_new, thresh = 0.60, n_fold = 5, iter = 1, deg = 3, DE_MEAN_RESP = T,
                                        pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = 5, tau = 0.75, spat_triang = kns_tri2,
                                        comp_mod = F, svfqm = T, lambda_seq = c(lambda[i],1))
  ########## Storing the results #############
  kns_svfqm_lam_QEP_tau0.25[2:6,i] = kns_qr_res_pnlty_tau0.25$svqfm$train_QEP
  kns_svfqm_lam_QEP_tau0.25[7:11,i] = kns_qr_res_pnlty_tau0.25$svqfm$test_QEP
  ############
  kns_svfqm_lam_QEP_tau0.5[2:6,i] = kns_qr_res_pnlty_tau0.5$svqfm$train_QEP
  kns_svfqm_lam_QEP_tau0.5[7:11,i] = kns_qr_res_pnlty_tau0.5$svqfm$test_QEP
  ########
  kns_svfqm_lam_QEP_tau0.75[2:6,i] = kns_qr_res_pnlty_tau0.75$svqfm$train_QEP
  kns_svfqm_lam_QEP_tau0.75[7:11,i] = kns_qr_res_pnlty_tau0.75$svqfm$test_QEP
}


# calculating the means
kns_svfqm_lambda_cv_tau0.25[,2] = apply(abs(kns_svfqm_lam_QEP_tau0.25[2:6,]),2,mean)
kns_svfqm_lambda_cv_tau0.25[,3] = apply(abs(kns_svfqm_lam_QEP_tau0.25[7:11,]),2,mean)
kns_svfqm_lambda_cv_tau0.5[,2] = apply(abs(kns_svfqm_lam_QEP_tau0.5[2:6,]),2,mean)
kns_svfqm_lambda_cv_tau0.5[,3] = apply(abs(kns_svfqm_lam_QEP_tau0.5[7:11,]),2,mean)
kns_svfqm_lambda_cv_tau0.75[,2] = apply(abs(kns_svfqm_lam_QEP_tau0.75[2:6,]),2,mean)
kns_svfqm_lambda_cv_tau0.75[,3] = apply(abs(kns_svfqm_lam_QEP_tau0.75[7:11,]),2,mean)

kns_svfqm_lambda_cv_tau0.25
kns_svfqm_lambda_cv_tau0.5
kns_svfqm_lambda_cv_tau0.75

# Saving the results for all folds
# save(file = kns_svfqm_lam_QEP_tau0.25, filename = "KnsSVFQM_PnltyParamCV_tau0.25.RData")

####### Creating a plot for the Lambda CV results #######
lambda =  c(0.00025,0.0005,0.001,0.0015,0.003 ,0.005,0.01,0.025,0.05,0.075,0.10,0.15,0.20)
### 25th quantile
kns_tau0.25_svfqm_lambda_df = data.frame(matrix(0, nrow = 2*length(lambda), ncol = 3))
kns_tau0.25_svfqm_lambda_df[,1] = lambda
kns_tau0.25_svfqm_lambda_df[1:13,2] = kns_svfqm_lambda_cv_tau0.25[,2] # train values
kns_tau0.25_svfqm_lambda_df[1:13,3] = "train"
kns_tau0.25_svfqm_lambda_df[14:26,2] = kns_svfqm_lambda_cv_tau0.25[,3] # test values
kns_tau0.25_svfqm_lambda_df[14:26,3] = "test"
colnames(kns_tau0.25_svfqm_lambda_df) = c("lambda","QEP","set")
### 50th quantile
kns_tau0.5_svfqm_lambda_df = data.frame(matrix(0, nrow = 2*length(lambda), ncol = 3))
kns_tau0.5_svfqm_lambda_df[,1] = lambda
kns_tau0.5_svfqm_lambda_df[1:13,2] = kns_svfqm_lambda_cv_tau0.5[,2] # train values
kns_tau0.5_svfqm_lambda_df[1:13,3] = "train"
kns_tau0.5_svfqm_lambda_df[14:26,2] = kns_svfqm_lambda_cv_tau0.5[,3] # test values
kns_tau0.5_svfqm_lambda_df[14:26,3] = "test"
colnames(kns_tau0.5_svfqm_lambda_df) = c("lambda","QEP","set")
### 75th quantile
kns_tau0.75_svfqm_lambda_df = data.frame(matrix(0, nrow = 2*length(lambda), ncol = 3))
kns_tau0.75_svfqm_lambda_df[,1] = lambda
kns_tau0.75_svfqm_lambda_df[1:13,2] = kns_svfqm_lambda_cv_tau0.75[,2] # train values
kns_tau0.75_svfqm_lambda_df[1:13,3] = "train"
kns_tau0.75_svfqm_lambda_df[14:26,2] = kns_svfqm_lambda_cv_tau0.75[,3] # test values
kns_tau0.75_svfqm_lambda_df[14:26,3] = "test"
colnames(kns_tau0.75_svfqm_lambda_df) = c("lambda","QEP","set")

ggplot(data = kns_tau0.25_svfqm_lambda_df, mapping = aes(x = lambda, y = QEP, color = set, linetype = set)) + geom_point() + 
  geom_line() + xlab("Nu") + ylab("QEP")

ggplot(data = kns_tau0.5_svfqm_lambda_df, mapping = aes(x = lambda, y = QEP, color = set, linetype = set)) + geom_point() + 
  geom_line() + xlab("Nu") + ylab("QEP")

ggplot(data = kns_tau0.75_svfqm_lambda_df, mapping = aes(x = lambda, y = QEP, color = set, linetype = set)) + geom_point() + 
  geom_line() + xlab("Nu") + ylab("QEP")


################################################################################################################################################
########################################## KANSAS DATA CROSS-VALIDATION FOR FPCA ############################################################### 
################################################################################################################################################
NUM_ITERS = 9
fld_lim = 5
# Creating lists for train and test mspe for both SVFQM and FQR
ag_data_svfqm_QEP_res = list()
ag_data_fqr_QEP_res = list()
num_scores_used_svfqm = list()
num_scores_used_fqr = list()
# perc_var can be adjusted for fine-tuning
perc_var = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
tau = c(0.25,0.5,0.75)
# Optimal penalty values 
lambda_vals = c(0.0005,0.0015,0.05)
tau_QEP_res_svfqm = list()
tau_QEP_res_fqr = list()
for(i in 1:length(tau)){
  lambda = lambda_vals[i] # the lambda corresponding to the tau'th quantile
  cat("The current quantile is ", tau[i], "\n") 
  for(j in 1:length(perc_var)){
    cat("  The current percentage of variance is ", perc_var[j], "\n") 
    # Organizing the matrices to store the train and test QEP for each percentage
    ag_data_QEP_svfqm = data.frame(matrix(0, nr = 45, nc = 2))
    colnames(ag_data_QEP_svfqm) = c("train QEP", "test QEP")
    ag_data_QEP_fqr = data.frame(matrix(0, nr = 45, nc = 2))
    colnames(ag_data_QEP_fqr) = c("train QEP", "test QEP")
    
    for (k in 1:NUM_ITERS) {
      ag.data.cv_otpt = ag.data_qr(fundat_all = kns_fd_new, nonfd = kns_nonfd_new, thresh = perc_var[j], n_fold = 5, iter = k, deg = 3, DE_MEAN_RESP = T,
                                   pred_vars = c("avgPRCP","irrig_prop"), reval = c(0,1), fld_lim = fld_lim, tau = tau[i], spat_triang = kns_tri2,
                                   comp_mod = T, svfqm = T, lambda_seq = c(lambda,1.35))
      # Updating SVFQM results
      ag_data_QEP_svfqm[((k-1)*5 + 1):(5*k),1] = ag.data.cv_otpt$svqfm$train_QEP
      ag_data_QEP_svfqm[((k-1)*5 + 1):(5*k),2] = ag.data.cv_otpt$svqfm$test_QEP
      # Updating FQR results
      ag_data_QEP_fqr[((k-1)*5 + 1):(5*k),1] = ag.data.cv_otpt$fqr$train_QEP
      ag_data_QEP_fqr[((k-1)*5 + 1):(5*k),2] = ag.data.cv_otpt$fqr$test_QEP
    }
    ag_data_svfqm_QEP_res[[j]] = ag_data_QEP_svfqm
    ag_data_fqr_QEP_res[[j]] = ag_data_QEP_fqr
    # The same number of scores is used for each repition so I only need to keep
    # track of it at the end
    num_scores_used_svfqm[[j]] = ag.data.cv_otpt$svqfm$num_used_harm
    num_scores_used_fqr[[j]] = ag.data.cv_otpt$fqr$num_used_harm
  }
  # Store results for the tau'th quantile
  tau_QEP_res_svfqm[[i]] = ag_data_svfqm_QEP_res
  tau_QEP_res_fqr[[i]] = ag_data_fqr_QEP_res
}

# saving the results 
###### Saving updated FQR CV results#######
#date = format(Sys.time(), "%m-%d-%Y")
# file_name = paste("UpdatedKnsFQR_CV_W_Pnlty_Results_PercVar0.15to0.90_0.95_tau_0.25_0.75", date, ".RData", sep = "")
# save(tau_QEP_res_fqr, file = file_name)
## Fine-tuning for updated FQR
# date = format(Sys.time(), "%m-%d-%Y")
# file_name = paste("UpdatedKnsFQR_FineTune_W_Pnlty_Results_tau_0.25_0.75", date, ".RData", sep = "")
# save(tau_QEP_res_fqr, file = file_name)# 
#### Saving updated SVFQM results with new lambda #######
# date = format(Sys.time(), "%m-%d-%Y")
# file_name = paste("UpdtLambdaKnsSVFQM_CV_W_Pnlty_Res_PercVar0.15to0.90_0.95_tau_0.25_0.50_0.75", date, ".RData", sep = "")
# save(tau_QEP_res_svfqm, file = file_name)

####### Calculating the average SEP for the percentiles and quantiles above ########
# perc_var can be adjusted for fine-tuning
perc_var = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
tau_0.25_avg_sep_svfqm = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.25_avg_sep_svfqm[,1] = perc_var
colnames(tau_0.25_avg_sep_svfqm) = c("perc_var","train_QEP", "test_QEP")
tau_0.25_avg_sep_fqr = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.25_avg_sep_fqr[,1] = perc_var
colnames(tau_0.25_avg_sep_fqr) = c("perc_var","train_QEP", "test_QEP")
tau_0.5_avg_sep_svfqm = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.5_avg_sep_svfqm[,1] = perc_var
colnames(tau_0.5_avg_sep_svfqm) = c("perc_var","train_QEP", "test_QEP")
tau_0.5_avg_sep_fqr = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.5_avg_sep_fqr[,1] = perc_var
colnames(tau_0.5_avg_sep_fqr) = c("perc_var","train_QEP", "test_QEP")
tau_0.75_avg_sep_svfqm = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.75_avg_sep_svfqm[,1] = perc_var
colnames(tau_0.75_avg_sep_svfqm) = c("perc_var","train_QEP", "test_QEP")
tau_0.75_avg_sep_fqr = matrix(0, nrow = length(perc_var), ncol = 3)
tau_0.75_avg_sep_fqr[,1] = perc_var
colnames(tau_0.75_avg_sep_fqr) = c("perc_var","train_QEP", "test_QEP")

# Calculating the averages
for(j in 1:length(perc_var)){
  tau_0.25_avg_sep_svfqm[j,2] = mean(abs(tau_QEP_res_svfqm[[1]][[j]]$`train QEP`))
  tau_0.25_avg_sep_svfqm[j,3] = mean(abs(tau_QEP_res_svfqm[[1]][[j]]$`test QEP`))
  tau_0.25_avg_sep_fqr[j,2] = mean(abs(tau_QEP_res_fqr[[1]][[j]]$`train QEP`))
  tau_0.25_avg_sep_fqr[j,3] = mean(abs(tau_QEP_res_fqr[[1]][[j]]$`test QEP`))

  tau_0.5_avg_sep_svfqm[j,2] = mean(abs(tau_QEP_res_svfqm[[2]][[j]]$`train QEP`))
  tau_0.5_avg_sep_svfqm[j,3] = mean(abs(tau_QEP_res_svfqm[[2]][[j]]$`test QEP`))
  tau_0.5_avg_sep_fqr[j,2] = mean(abs(tau_QEP_res_fqr[[2]][[j]]$`train QEP`))
  tau_0.5_avg_sep_fqr[j,3] = mean(abs(tau_QEP_res_fqr[[2]][[j]]$`test QEP`))

  tau_0.75_avg_sep_svfqm[j,2] = mean(abs(tau_QEP_res_svfqm[[3]][[j]]$`train QEP`))
  tau_0.75_avg_sep_svfqm[j,3] = mean(abs(tau_QEP_res_svfqm[[3]][[j]]$`test QEP`))
  tau_0.75_avg_sep_fqr[j,2] = mean(abs(tau_QEP_res_fqr[[3]][[j]]$`train QEP`))
  tau_0.75_avg_sep_fqr[j,3] = mean(abs(tau_QEP_res_fqr[[3]][[j]]$`test QEP`))
}

tau_0.25_avg_sep_svfqm
tau_0.25_avg_sep_fqr

tau_0.5_avg_sep_svfqm
tau_0.5_avg_sep_fqr

tau_0.75_avg_sep_svfqm
tau_0.75_avg_sep_fqr

num_scores_used_svfqm
num_scores_used_fqr

#################################### Kansas SVFQM and FQR CV Perc Var Plots ####################################
library(tidyverse)
###### Results for tau = 0.25 ######
##### SVFQM results
# Setting it up for SVFQM with the first couple rows correspond to train and the last couple rows for test. 
# This is the way I set it up in the dataframe anyways which is easier to work with for ggplot
kns_QEP_scores_cv_svfqm_tau0.25 = matrix(0, nrow = 14, 2)
kns_QEP_scores_cv_svfqm_tau0.25[1:7,1] = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
kns_QEP_scores_cv_svfqm_tau0.25[8:14,1] = kns_QEP_scores_cv_svfqm_tau0.25[1:7,1]
kns_QEP_scores_cv_svfqm_tau0.25[1:7,2] = c(0.0039,0.004,0.0036,0.0046,0.0036,0.0029,0.0040) # the train QEP
kns_QEP_scores_cv_svfqm_tau0.25[8:14,2] = c(0.0228,0.0234,0.0223,0.0262,0.0214,0.0236,0.0215)  # test QEP
kns_QEP_scores_cv_svfqm_tau0.25_df = data.frame(kns_QEP_scores_cv_svfqm_tau0.25)
kns_QEP_scores_cv_svfqm_tau0.25_df[1:7,"Set"] = "train"
kns_QEP_scores_cv_svfqm_tau0.25_df[8:14,"Set"] = "test"
colnames(kns_QEP_scores_cv_svfqm_tau0.25_df) = c("PercVar","QEP", "Set")
ggplot(data = kns_QEP_scores_cv_svfqm_tau0.25_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")

##### FQR results 
kns_QEP_scores_cv_fqr_tau0.25 = matrix(0, nrow = 14, 2)
kns_QEP_scores_cv_fqr_tau0.25[1:7,1] = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
kns_QEP_scores_cv_fqr_tau0.25[8:14,1] = kns_QEP_scores_cv_fqr_tau0.25[1:7,1]
kns_QEP_scores_cv_fqr_tau0.25[1:7,2] = c(0.0015,0.0017,0.0025,0.0032,0.0037,0.0047,0.0072) # the train QEP
kns_QEP_scores_cv_fqr_tau0.25[8:14,2] = c(0.0223,0.0223,0.0244,0.0207,0.0225,0.0271,0.023)  # test QEP
kns_QEP_scores_cv_fqr_tau0.25_df = data.frame(kns_QEP_scores_cv_fqr_tau0.25)
kns_QEP_scores_cv_fqr_tau0.25_df[1:7,"Set"] = "train"
kns_QEP_scores_cv_fqr_tau0.25_df[8:14,"Set"] = "test"
colnames(kns_QEP_scores_cv_fqr_tau0.25_df) = c("PercVar","QEP", "Set")
ggplot(data = kns_QEP_scores_cv_fqr_tau0.25_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")


###### Results for tau = 0.50 ######

##### SVFQM results
kns_QEP_scores_cv_svfqm_tau0.50 = matrix(0, nrow = 14, 2)
kns_QEP_scores_cv_svfqm_tau0.50[1:7,1] = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
kns_QEP_scores_cv_svfqm_tau0.50[8:14,1] = kns_QEP_scores_cv_svfqm_tau0.50[1:7,1]
kns_QEP_scores_cv_svfqm_tau0.50[1:7,2] = c(0.0082,0.0081,0.0086,0.0077,0.0078,0.0090,0.0088) # the train QEP
kns_QEP_scores_cv_svfqm_tau0.50[8:14,2] = c(0.0275,0.0260,0.0249,0.0263,0.0232,0.0264,0.0262)  # test QEP
kns_QEP_scores_cv_svfqm_tau0.50_df = data.frame(kns_QEP_scores_cv_svfqm_tau0.50)
kns_QEP_scores_cv_svfqm_tau0.50_df[1:7,"Set"] = "train"
kns_QEP_scores_cv_svfqm_tau0.50_df[8:14,"Set"] = "test"
colnames(kns_QEP_scores_cv_svfqm_tau0.50_df) = c("PercVar","QEP", "Set")
ggplot(data = kns_QEP_scores_cv_svfqm_tau0.50_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")

##### FQR results 
kns_QEP_scores_cv_fqr_tau0.50 = matrix(0, nrow = 18, 2)
kns_QEP_scores_cv_fqr_tau0.50[1:9,1] = c(0.15,0.30,0.45,0.60,0.75,0.80,0.85,0.90,0.95)
kns_QEP_scores_cv_fqr_tau0.50[10:18,1] = kns_QEP_scores_cv_fqr_tau0.50[1:9,1]
kns_QEP_scores_cv_fqr_tau0.50[1:9,2] = c(0.0014,0.0020,0.0022,0.0032,0.0043,0.0046,0.0047,0.0051,0.0082) # the train QEP
kns_QEP_scores_cv_fqr_tau0.50[10:18,2] = c(0.0248,0.0239,0.0237,0.0242,0.0233,0.0238,0.0249,0.0195,0.0204) # test QEP
kns_QEP_scores_cv_fqr_tau0.50_df = data.frame(kns_QEP_scores_cv_fqr_tau0.50)
kns_QEP_scores_cv_fqr_tau0.50_df[1:9,"Set"] = "train"
kns_QEP_scores_cv_fqr_tau0.50_df[10:18,"Set"] = "test"
colnames(kns_QEP_scores_cv_fqr_tau0.50_df) = c("PercVar","QEP", "Set")
ggplot(data = kns_QEP_scores_cv_fqr_tau0.50_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")


###### Results for tau = 0.75 ######

#### SVFQM results
kns_QEP_scores_cv_svfqm_tau0.75 = matrix(0, nrow = 14, 2)
kns_QEP_scores_cv_svfqm_tau0.75[1:7,1] = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
kns_QEP_scores_cv_svfqm_tau0.75[8:14,1] = kns_QEP_scores_cv_svfqm_tau0.75[1:7,1]
kns_QEP_scores_cv_svfqm_tau0.75[1:7,2] = c(0.0069,0.0049,0.0046,0.0044,0.0059,0.0071,0.0063) # the train QEP
kns_QEP_scores_cv_svfqm_tau0.75[8:14,2] = c(0.0209,0.0207,0.0211,0.0209,0.0219,0.0207,0.0199)  # test QEP
kns_QEP_scores_cv_svfqm_tau0.75_df = data.frame(kns_QEP_scores_cv_svfqm_tau0.75)
kns_QEP_scores_cv_svfqm_tau0.75_df[1:7,"Set"] = "train"
kns_QEP_scores_cv_svfqm_tau0.75_df[8:14,"Set"] = "test"
colnames(kns_QEP_scores_cv_svfqm_tau0.75_df) = c("PercVar","QEP", "Set")
ggplot(data = kns_QEP_scores_cv_svfqm_tau0.75_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")

##### FQR results 
kns_QEP_scores_cv_fqr_tau0.75 = matrix(0, nrow = 14, 2)
kns_QEP_scores_cv_fqr_tau0.75[1:7,1] = c(0.15,0.30,0.45,0.60,0.75,0.90,0.95)
kns_QEP_scores_cv_fqr_tau0.75[8:14,1] = kns_QEP_scores_cv_fqr_tau0.75[1:7,1]
kns_QEP_scores_cv_fqr_tau0.75[1:7,2] = c(0.0012,0.0018,0.0023,0.0031,0.0037,0.0051,0.0065) # the train QEP
kns_QEP_scores_cv_fqr_tau0.75[8:14,2] = c(0.0181,0.0185,0.0185,0.0171,0.0179,0.0190,0.0206)  # test QEP
kns_QEP_scores_cv_fqr_tau0.75_df = data.frame(kns_QEP_scores_cv_fqr_tau0.75)
kns_QEP_scores_cv_fqr_tau0.75_df[1:7,"Set"] = "train"
kns_QEP_scores_cv_fqr_tau0.75_df[8:14,"Set"] = "test"
colnames(kns_QEP_scores_cv_fqr_tau0.75_df) = c("PercVar","QEP", "Set")
ggplot(data = kns_QEP_scores_cv_fqr_tau0.75_df, mapping = aes(x = PercVar, y = QEP, color = Set, linetype = Set)) + geom_point() + 
  geom_line() + xlab("Variance Percentage") + ylab("QEP")



















