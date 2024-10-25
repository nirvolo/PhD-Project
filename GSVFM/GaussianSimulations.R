##############################################################################################################################################################s
# Simulation Trials
##############################################################################################################################################################
# Assume opening the project at root directory
source("GSVFM/sim_funcs/GSVFM_SimSettings.R")

########## Tuning Parameters ################
N = 100 # number of time points
fd_err = sqrt(0.2) # sd of error (functional predictor)
# Matern corr. fxn parameter that controls rate of decay. Dividing by 111 km since
# 50 is in kilometer units, not lat and long units. The distance that's calculated
# in the simulation is between lat and long units.
phi = 50/111
yrs = 10 
flds = 5 
n_iter = 25 # number of iterations

#### **Example of running the simulation function for 1 iteration** #####
sim_ex_run = sim_fxn_new(N = N, phi = phi, v = fd_err, sp_locs = scl_locs, years = yrs, n_folds = flds, 
                         coef_comp = "complex", snr = 3.3, comp_model = T, num_iter = 1, spat_corr = T, d = 3, tri_fin = 2,
                         reval = c(0,1), log.lam = c(-15.5,-8), models = c("svfm","plfam","flm"), family = gaussian())



##############################################################################################
######### Cross-validation for different combinations of degree and triangulation ###########
##############################################################################################
# There are 6 combinations of degree and triangulation. Using degree's = 2, 3 and 
# triangulation fineness n = 1.2, 2
smooth_params = data.frame(matrix(c(2,1.2,3,1.2,2,2,3,2), nrow = 4, ncol = 2, byrow = T))
colnames(smooth_params) = c("degree", "tri_fine")
smooth_params_otpt_lst_spat_corr = list()
for(j in 1:4){
  deg = smooth_params$degree[j]
  tri_fine = smooth_params$tri_fine[j]
  cat("The degree is", deg, "and triangulation smoothness is", tri_fine, "\n")
  # With spatial correlation
  smooth_params_otpt_lst_spat_corr[[j]] = sim_fxn_new(N = N, phi = phi, v = fd_err, sp_locs = scl_locs, years = yrs, n_folds = flds, 
                                            coef_comp = "complex", snr = 3.3, comp_model = F, num_iter = n_iter, spat_corr = T, d = deg, tri_fin = tri_fine,
                                            reval = c(0,1), log.lam = c(-15.5,-8), models = c("svfm"), family = gaussian())
}

# Saving the results for the simulation trials for (deg, num_tri)
#save(smooth_params_otpt_lst_spat_corr, file = "SmoothParamsTrials_SpatCorrFINALVERSION_June13_2024.RData")
#save(smooth_params_otpt_lst_no_spat_corr, file = "SmoothParamsTrials_NoSpatCorr.RData")

# Summarizing the results from the above simulations (with spatial correlation)
svfm_mspe_deg_two_three = matrix(NA, nr = 4, nc = 125)
svfm_mspe_var_deg_two_three = matrix(NA, nr = 4, nc = 125)
svfm_runtime_deg_two_three = matrix(NA, nr = 4, nc = 125)
avg_svfm_mise_deg_two_three = data.frame(matrix(NA, nr = 4, nc = 7))
for(i in 1:4){
  svfm_mspe_deg_two_three[i,] = smooth_params_otpt_lst_spat_corr[[i]]$svfm$svfm_mspe$test_MSE
  svfm_mspe_var_deg_two_three[i,] = smooth_params_otpt_lst_spat_corr[[i]]$svfm$svfm_mspe$test_mse_var
  svfm_runtime_deg_two_three[i,] = smooth_params_otpt_lst_spat_corr[[i]]$svfm$svfm_train_time
  avg_svfm_mise_deg_two_three[i,] = apply(smooth_params_otpt_lst_spat_corr[[i]]$svfm$beta_mise, 2 , mean)
}

# Calculating the mean for each setting
avg_svfm_mspe_deg_two_three = apply(svfm_mspe_deg_two_three, 1, mean)
avg_svfm_mspe_var_deg_two_three = apply(svfm_mspe_var_deg_two_three, 1, mean)
avg_svfm_runtime_deg_two_three = apply(svfm_runtime_deg_two_three, 1, mean)
avg_svfm_mspe_deg_two_three
avg_svfm_mspe_var_deg_two_three
avg_svfm_runtime_deg_two_three
avg_svfm_mise_deg_two_three


############################################################################################################################################
# Running all simulation settings for GSVFM, PLFAM, and FLM
############################################################################################################################################
sim_settings = matrix(c(F, "const", F, "basic", F, "complex", T, "const", T, "basic", T, "complex"), nrow = 6, ncol = 2, byrow = T)
# A list to store the output from each simulation trial
sim_otpt_lst_settings = list(set1 = NULL, set2 = NULL, set3 = NULL, set4 = NULL, set5 = NULL, set6 = NULL)
# Iterating over the different settings
yrs = 10 # Change to 10
flds = 5 # Change to 5
n_iter = 25 # number of iterations
for(i in 1:6){
  crnt_set = sim_settings[i,]
  cat("Setting", i, "\n")
  sim_otpt_lst_settings[[i]] = sim_fxn_new(N = N, phi = phi, v = fd_err, sp_locs = scl_locs, years = yrs, n_folds = flds, 
                                      coef_comp = crnt_set[2], snr = 3.3, comp_model = T, num_iter = n_iter, spat_corr = crnt_set[1], d = 3, tri_fin = 2,
                                      reval = c(0,1), log.lam = c(-15.5,-8), models = c("svfm","plfam","flm"), family = gaussian())
}
#save(sim_otpt_lst_settings, file = "SimOtptResultsAllSettingsMarch31_2024.RData")
#save(sim_otpt_lst_settings, file = "SimOtptResultsAllSettingsFINAL_VERSION_Feb1_2024.RData")

sim_otpt_lst_settings[[1]]$svfm$svfm_train_time
all_fxns_svfm_mspe = matrix(NA, nr = 6, nc = 5)
all_fxns_svfm_runtime = matrix(NA, nr = 6, nc = 1)
all_fxns_plfam_mspe = matrix(NA, nr = 6, nc = 5)
all_fxns_plfam_runtime = matrix(NA, nr = 6, nc = 1)
all_fxns_flm_mspe = matrix(NA, nr = 6, nc = 5)
all_fxns_flm_runtime = matrix(NA, nr = 6, nc = 1)
all_fxns_svfm_mise = matrix(NA, nr = 6, nc = 7)
for(i in 1:6){
  all_fxns_svfm_mspe[i,] = apply(sim_otpt_lst_settings[[i]]$svfm$svfm_mspe, 2, mean)
  all_fxns_svfm_runtime[i,] = mean(sim_otpt_lst_settings[[i]]$svfm$svfm_train_time)
  all_fxns_plfam_mspe[i,] = apply(sim_otpt_lst_settings[[i]]$plfam$plfam_mspe, 2, mean)
  all_fxns_plfam_runtime[i,] = mean(sim_otpt_lst_settings[[i]]$plfam$plfam_train_time)
  all_fxns_svfm_mise[i,] = apply(sim_otpt_lst_settings[[i]]$svfm$beta_mise, 2, mean)
  all_fxns_flm_mspe[i,] = apply(sim_otpt_lst_settings[[i]]$flm$flm_mspe, 2, mean)
  all_fxns_flm_runtime[i,] = mean(sim_otpt_lst_settings[[i]]$flm$flm_train_time)
}

colnames(all_fxns_svfm_mspe) = colnames(sim_otpt_lst_settings[[1]]$svfm$svfm_mspe)
colnames(all_fxns_plfam_mspe) = colnames(sim_otpt_lst_settings[[1]]$svfm$svfm_mspe)
colnames(all_fxns_flm_mspe) = colnames(sim_otpt_lst_settings[[1]]$svfm$svfm_mspe)
all_fxns_svfm_mspe
all_fxns_svfm_runtime
all_fxns_plfam_mspe
all_fxns_plfam_runtime
all_fxns_flm_mspe
all_fxns_flm_runtime
all_fxns_svfm_mise

###################### Boxplots of for the MISE for the complex setting and spatial correlation #####################
sim_mspe_df_comp = data.frame(matrix(0,nrow = 375, ncol = 2))
sim_mspe_df_comp[1:125,1] = "SVFM"
sim_mspe_df_comp[1:125,2] = sim_otpt_lst_settings$set6$svfm$svfm_mspe$test_MSE
sim_mspe_df_comp[126:250,1] = "PLFAM"
sim_mspe_df_comp[126:250,2] = sim_otpt_lst_settings$set6$plfam$plfam_mspe$test_MSE
sim_mspe_df_comp[251:375,1] = "FLM"
sim_mspe_df_comp[251:375,2] = sim_otpt_lst_settings$set6$flm$flm_mspe$test_MSE
colnames(sim_mspe_df_comp) = c("Model","MSPE")
sim_mspe_df_comp$Model = factor(sim_mspe_df_comp$Model, levels = c("SVFM","PLFAM","FLM"), ordered = TRUE) d
ggplot(data = sim_mspe_df_comp, aes(x = sim_mspe_df_comp$Model, y = sim_mspe_df_comp$MSPE)) +
  geom_boxplot() + 
  xlab("Model") +
  ylab("MSPE") +
  ggtitle("Model Comparison (Complex and SC)")
  









