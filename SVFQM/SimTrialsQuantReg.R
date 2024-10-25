source("SVFQM/sim_funcs/SimFxnQuantReg.R")

##### Setting up parameters for simulation #####
N = 100 # number of time points
fd_err = sqrt(0.2) # sd of error (functional predictor)
# Matern corr. fxn parameter that controls rate of decay. Dividing by 111 km since
# 50 is in kilometer units, not lat and long units. The distance that's calculated
# in the simulation is between lat and long units.
phi = 50/111
# Iterating over the different settings
yrs = 10 # Change to 10
flds = 4 # Change to 5
n_iter = 25 # number of iterations

##### Example of running one iteration of the simulation ######
sim_qr_ex = sim_fxn_qr(N = N, phi = phi, v = fd_err, sp_locs = scl_locs, years = 10, n_folds = 5, 
                         coef_comp = "complex", snr = 3.3, num_iter = 1, spat_corr = F, d = 3, tri_fin = 2,
                         reval = c(0,1), log.lam = c(-15.5,-8), models = c("svqfm","fqr"), tau = 0.5, comp_model = T, scl = T)

###### Running 25 repetitions for all settings using different values of tau #########
sim_settings = matrix(c(F, "const", F, "basic", F, "complex", T, "const", T, "basic", T, "complex"), nrow = 6, ncol = 2, byrow = T)
quant_vals = c(0.25,0.5,0.75)
# List to store the results for each value of tau
qr_sim_pred_err = list("tau_0.25","tau_0.5","tau_0.75")

for(i in 1:length(quant_vals)){
  tau = quant_vals[i]
  cat("The value of tau is", tau,"\n")
  sim_otpt_lst_settings = list(set1 = NULL, set2 = NULL, set3 = NULL, set4 = NULL, set5 = NULL, set6 = NULL)
  for(j in 1:nrow(sim_settings)){
    crnt_set = sim_settings[j,]
    sc = crnt_set[1]
    fxn = crnt_set[2]
    if(sc == T) cat("The current setting is spatial correlation and", fxn, "functions", "\n")
    if(sc == F) cat("The current setting is no spatial correlation and", fxn, "functions", "\n")
    sim_qr_res = sim_fxn_qr(N = N, phi = phi, v = fd_err, sp_locs = scl_locs, years = yrs, n_folds = flds, 
                           coef_comp = fxn, snr = 3.3, num_iter = n_iter, spat_corr = sc, d = 3, tri_fin = 2,
                           reval = c(0,1), log.lam = c(-15.5,-8), models = c("svqfm","fqr"), tau = tau, comp_model = T, scl = T)
    sim_otpt_lst_settings[[j]] = sim_qr_res
  }
  # The results for all settings for the ith quantile
  qr_sim_pred_err[[i]] = sim_otpt_lst_settings
}

#save(qr_sim_pred_err, file = "QuantileSimResultsAllSettingsWithPenalty_tau_0.25_0.5_0.75_August20_2024.RData")


### Calculating the train and test pred error averages using the data above
# Creating a matrix for each value of tau with each row corresponding to the setting
# and each column corresponding to train and test respectively
sim_res_tau0.25_svqfm = matrix(0, nrow = 6, ncol = 2)
colnames(sim_res_tau0.25_svqfm) = c("train error", "test error")
sim_res_tau0.50_svqfm = matrix(0, nrow = 6, ncol = 2)
colnames(sim_res_tau0.50_svqfm) = c("train error", "test error")
sim_res_tau0.75_svqfm = matrix(0, nrow = 6, ncol = 2)
colnames(sim_res_tau0.75_svqfm) = c("train error", "test error")
sim_res_tau0.25_fqr = matrix(0, nrow = 6, ncol = 2)
colnames(sim_res_tau0.25_fqr) = c("train error", "test error")
sim_res_tau0.50_fqr = matrix(0, nrow = 6, ncol = 2)
colnames(sim_res_tau0.50_fqr) = c("train error", "test error")
sim_res_tau0.75_fqr = matrix(0, nrow = 6, ncol = 2)
colnames(sim_res_tau0.75_fqr) = c("train error", "test error")

for(j in 1:6){
  sim_res_tau0.25_svqfm[j,1] = mean(qr_sim_pred_err[[1]][[j]]$svqfm$train_mse)
  sim_res_tau0.25_svqfm[j,2] = mean(qr_sim_pred_err[[1]][[j]]$svqfm$test_mse)
  sim_res_tau0.50_svqfm[j,1] = mean(qr_sim_pred_err[[2]][[j]]$svqfm$train_mse)
  sim_res_tau0.50_svqfm[j,2] = mean(qr_sim_pred_err[[2]][[j]]$svqfm$test_mse)
  sim_res_tau0.75_svqfm[j,1] = mean(qr_sim_pred_err[[3]][[j]]$svqfm$train_mse)
  sim_res_tau0.75_svqfm[j,2] = mean(qr_sim_pred_err[[3]][[j]]$svqfm$test_mse)
}

sim_res_tau0.25_svqfm = cbind(sim_res_tau0.25_svqfm, sim_settings)
sim_res_tau0.25_fqr = cbind(sim_res_tau0.25_fqr, sim_settings)

sim_res_tau0.50_svqfm = cbind(sim_res_tau0.50_svqfm, sim_settings)
sim_res_tau0.50_fqr = cbind(sim_res_tau0.50_fqr, sim_settings)

sim_res_tau0.75_svqfm = cbind(sim_res_tau0.75_svqfm, sim_settings)
sim_res_tau0.75_fqr = cbind(sim_res_tau0.75_fqr, sim_settings)

sim_res_tau0.25_svqfm
sim_res_tau0.25_fqr

sim_res_tau0.50_svqfm
sim_res_tau0.50_fqr

sim_res_tau0.75_svqfm
sim_res_tau0.75_fqr














