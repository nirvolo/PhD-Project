source("GSVFM/sim_funcs/GSVFM_SimSettings.R")
N = 100 # number of time points
fd_err = sqrt(0.2) # sd of error (functional predictor)
# Matern corr. fxn parameter that controls rate of decay. Dividing by 111 km since
# 50 is in kilometer units, not lat and long units. The distance that's calculated
# in the simulation is between lat and long units.
phi = 50/111
# Iterating over the different settings
yrs = 10
flds = 5 
SNR = 3.3
#### **Example running the simulation with poisson response** #####
sim_basic_pois_resp = sim_fxn_new(N = N, phi = phi, v = fd_err, sp_locs = scl_locs, years = yrs, n_folds = flds, fld_lim = 1,
                                 coef_comp = "basic", snr = SNR, comp_model = T, num_iter = 1, spat_corr = T, d = 3, tri_fin = 2, models = c("svfm","gflm"), 
                                 reval = c(0,1), log.lam = c(-15.5,-8), family = poisson())

################ Running the Poisson simulations for all settings for GSVFM ##################
sim_settings = matrix(c(F, "const", F, "basic", F, "complex", T, "const", T, "basic", T, "complex"), nrow = 6, ncol = 2, byrow = T)
# A list to store the output from each simulation trial
sim_otpt_lst_settings = list(set1 = NULL, set2 = NULL, set3 = NULL, set4 = NULL, set5 = NULL, set6 = NULL)
# Iterating over the different settings
yrs = 10 # Change to 10
flds = 5 # Change to 5
n_iter = 25 # number of iterations
N = 100 # number of time points
fd_err = sqrt(0.2) # sd of error (functional predictor)
# Matern corr. fxn parameter that controls rate of decay. Dividing by 111 km since
# 50 is in kilometer units, not lat and long units. The distance that's calculated
# in the simulation is between lat and long units.
phi = 50/111
# Iterating over the different settings
yrs = 10 #
flds = 5
SNR = 3.3
for(i in 1:6){
  crnt_set = sim_settings[i,]
  cat("Setting", i, "\n")
  sim_otpt_lst_settings[[i]] = sim_fxn_new(N = N, phi = phi, v = fd_err, sp_locs = scl_locs, years = yrs, n_folds = flds, fld_lim = 5,
                                           coef_comp = crnt_set[2], snr = SNR, comp_model = T, num_iter = n_iter, spat_corr = crnt_set[1], d = 3, tri_fin = 2, models = c("svfm"), 
                                           reval = c(0,1), log.lam = c(-15.5,-8), family = poisson())
}
# Saving the file with the current date
date = format(Sys.time(), "%m-%d-%Y")
file_name_pois_sim = paste("PoissonResponse_SimOtptResultsAllSettings_FINAL_VERSION", date, ".RData", sep = "")
#save(sim_otpt_lst_settings, file = file_name_pois_sim)
# Calculating the average errors
poisson_all_fxns_svfm_mspe = matrix(NA, nr = 6, nc = 5)
poisson_all_fxns_svfm_runtime = matrix(NA, nr = 6, nc = 1)
poisson_all_fxns_svfm_mise = matrix(NA, nr = 6, nc = 7)
for(i in 1:6){
  poisson_all_fxns_svfm_mspe[i,] = apply(sim_otpt_lst_settings[[i]]$svfm$svfm_mspe, 2, mean)
  poisson_all_fxns_svfm_runtime[i,] = mean(sim_otpt_lst_settings[[i]]$svfm$svfm_train_time)
  poisson_all_fxns_svfm_mise[i,] = apply(sim_otpt_lst_settings[[i]]$svfm$beta_mise, 2, mean)
}
colnames(poisson_all_fxns_svfm_mspe) = colnames(sim_otpt_lst_settings[[1]]$svfm$svfm_mspe)
poisson_all_fxns_svfm_mspe


#################################### Poisson simulations for Generalized FLM #################################### 
######## Running all simulation settings for the Poisson FLM #########
sim_settings = matrix(c(F, "const", F, "basic", F, "complex", T, "const", T, "basic", T, "complex"), nrow = 6, ncol = 2, byrow = T)
sim_otpt_lst_settings_pois_flm = list(set1 = NULL, set2 = NULL, set3 = NULL, set4 = NULL, set5 = NULL, set6 = NULL)
# Iterating over the different settings
yrs = 10 # Change to 10
flds = 5 # Change to 5
n_iter = 25 # number of iterations
N = 100 # number of time points
fd_err = sqrt(0.2) # sd of error (functional predictor)
# Matern corr. fxn parameter that controls rate of decay. Dividing by 111 km since
# 50 is in kilometer units, not lat and long units. The distance that's calculated
# in the simulation is between lat and long units.
phi = 50/111
# Iterating over the different settings
yrs = 10 # Change to 10
flds = 5 # Change to 5
SNR = 3.3
start = Sys.time()
for(i in 1:6){
  crnt_set = sim_settings[i,]
  cat("Setting", i, "\n")
  sim_otpt_lst_settings_pois_flm[[i]] = sim_fxn_new(N = N, phi = phi, v = fd_err, sp_locs = scl_locs, years = yrs, n_folds = flds, fld_lim = 5,
                                           coef_comp = crnt_set[2], snr = SNR, comp_model = T, num_iter = n_iter, spat_corr = crnt_set[1], d = 3, tri_fin = 2, models = c("gflm"), 
                                           reval = c(0,1), log.lam = c(-15.5,-8), family = poisson())
}
# Saving the output
date = format(Sys.time(), "%m-%d-%Y")
file_name_pois_sim_flm = paste("PoissonResponseFLM_SimOtptResultsAllSettings_FINAL_VERSION", date, ".RData", sep = "")
#save(sim_otpt_lst_settings_pois_flm, file = file_name_pois_sim_flm)
# Calculating the average error
sim_otpt_lst_settings_pois_flm[[1]]$gflm$gflm_mspe
poisson_all_fxns_flm_mspe = matrix(NA, nr = 6, nc = 5)
poisson_all_fxns_flm_runtime = matrix(NA, nr = 6, nc = 1)
poisson_all_fxns_flm_mise = matrix(NA, nr = 6, nc = 7)
for(i in 1:6){
  poisson_all_fxns_flm_mspe[i,] = apply(sim_otpt_lst_settings_pois_flm[[i]]$gflm$gflm_mspe, 2, mean)
  poisson_all_fxns_flm_runtime[i,] = mean(sim_otpt_lst_settings_pois_flm[[i]]$gflm$pois_flm_train_time)
  #poisson_all_fxns_flm_mise[i,] = apply(sim_otpt_lst_settings_pois_flm[[i]]$gflm$beta_mise, 2, mean)
}

colnames(poisson_all_fxns_flm_mspe) = colnames(sim_otpt_lst_settings_pois_flm[[1]]$gflm$gflm_mspe)
poisson_all_fxns_flm_mspe














