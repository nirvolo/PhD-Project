source("GSVFM/util_funcs/plfam_yehua.R")

run_model_comp = function(train.lst_fd, test.lst_fd, nonfd.train, nonfd.test, fpca_thresh = 0.95, models, num_harm = 5){
  ############ Converting the training and testing functional data from lists to matrices #################
  n_cnty = length(train.lst_fd) # number of counties in train and test
  n_train = nrow(nonfd.train)
  n_test = nrow(nonfd.test)
  #n_vars = ncol(nonfd.train)
  n_knots = nrow(train.lst_fd[[1]]$coefs[,,1])
  # Arrays for to store the coefficient matrices for train and test
  arry_train = array(dim = c(n_knots, n_train, 2))
  arry_test = array(dim = c(n_knots, n_test, 2))
  # Initializing the dataframes by adding the first counties coefficients
  x1_train = train.lst_fd[[1]]$coefs[,,1]
  x2_train = train.lst_fd[[1]]$coefs[,,2]
  x1_test = test.lst_fd[[1]]$coefs[,,1]
  x2_test = test.lst_fd[[1]]$coefs[,,2]
  
  # Adding the coefficient matrices for each county for both variables
  for(i in 2:n_cnty){
    fd_cnty_train = train.lst_fd[[i]]$coefs
    fd_cnty_test = test.lst_fd[[i]]$coefs
    x1_train = cbind(x1_train, fd_cnty_train[,,1])
    x2_train = cbind(x2_train, fd_cnty_train[,,2])
    x1_test = cbind(x1_test, fd_cnty_test[,,1])
    x2_test = cbind(x2_test, fd_cnty_test[,,2])
  }
  
  # Adding the x1 and x2 matrices to the train and test arrays
  arry_train[,,1] = x1_train
  arry_train[,,2] = x2_train
  arry_test[,,1] = x1_test
  arry_test[,,2] = x2_test
  # Creating the combined train and test fd object for fpca
  basis_fxns = train.lst_fd[[1]]$basis
  fd.train = fd(arry_train, basis_fxns)
  fd.test = fd(arry_test, basis_fxns)
  
  ############ Fitting the models and prediction #####################
  # The response and linear predictor variable
  lin_var.train = nonfd.train[,"z"]
  lin_var.test = nonfd.test[,"z"]
  resp.train = nonfd.train[,"y"]
  resp.test = nonfd.test[,"y"]
  
  #####################
  # PLFAM with COSSO
  #####################
  # The standardized pc scores to be used in the PLFAM model
  # In the next line we use the default number of harmonics in get.zeta
  # which is based on 99.9% variation explained. But in my model, I use only
  # 5 harmonics.
  plfam_fpca_start = unclass(Sys.time()) # unclass to convert to double in seconds
  std.scores_train = get.zeta(fdobj = fd.train, nharm = num_harm)$zeta
  t_fpca = unclass(Sys.time()) - plfam_fpca_start
  std.scores_test = get.zeta2(nfobj = fd.test, fdobj = fd.train, nharm = num_harm)$zeta

  # List to store the output from the different models 
  comp_mod_mspe = list()
  if("plfam" %in% models){
    # The PLFAM model fit using the training data. Measuring the run time of the model
    plfam_fit_start = unclass(Sys.time())
    plfam_fit = plcosso(y = resp.train, U = lin_var.train, zeta = std.scores_train)
    t_cosso = unclass(Sys.time()) - plfam_fit_start
    plfam_run_time = t_cosso + t_fpca
    # Performing prediction using the plfam model for train and test
    plfam_pred_train = predict.plcosso(obj = plfam_fit, Unew = lin_var.train, zetanew = std.scores_train)
    plfam_pred_test = predict.plcosso(obj = plfam_fit, Unew = lin_var.test, zetanew = std.scores_test)
    
    # Calculating the different types of MSPE for the PLFAM
    plfam_train_mse = mean((resp.train-plfam_pred_train)^2)
    plfam_mspe = mean((resp.test-plfam_pred_test)^2)
    plfam_mse_perc = mean((resp.test - plfam_pred_test)^2/(resp.test)^2) # normalizing by the square y
    plfam_mse_var = mean((resp.test - plfam_pred_test)^2)/var(resp.test) # normalizing by the variance of y
    plfam_mse_mean_abs = mean((resp.test - plfam_pred_test)^2)/mean(abs(resp.test)) # normalizing by the mean of the absolute value of y
    plfam_mse_mean_sqr = mean((resp.test - plfam_pred_test)^2)/mean(resp.test^2) # normalizing by the mean of the absolute value of y
    plfam_mspe_lst = list("plfam_mspe" = plfam_mspe, "plfam_mse_perc" = plfam_mse_perc, "plfam_mse_var" = plfam_mse_var,
                          "plfam_mean_abs" = plfam_mse_mean_abs, "plfam_mse_mean_sqr" = plfam_mse_mean_sqr, "plfam_run_time" = plfam_run_time, 
                          "plfam_train_MSE" = plfam_train_mse)
    comp_mod_mspe[["plfam"]] = plfam_mspe_lst
  }
  
  ##### FLM ################
  if("flm" %in% models){
    flm_fit_start = unclass(Sys.time())
    flm_fit <- myridge(y = resp.train, U = lin_var.train, R = std.scores_train)
    t_flm = unclass(Sys.time()) - flm_fit_start
    flm_run_time = t_flm + t_fpca
    flm_pred_train = as.vector(cbind(1, lin_var.train, std.scores_train) %*% flm_fit$beta) 
    flm_train_mse = mean((resp.train-flm_pred_train)^2)
    flm_pred_test = as.vector(cbind(1, lin_var.test, std.scores_test) %*% flm_fit$beta)
    # Calculating the different types of MSPE for the FLM
    flm_mspe = mean((resp.test-flm_pred_test)^2)
    flm_mse_perc = mean((resp.test - flm_pred_test)^2/(resp.test)^2) # normalizing by the square y
    flm_mse_var = mean((resp.test - flm_pred_test)^2)/var(resp.test) # normalizing by the variance of y
    flm_mse_mean_abs = mean((resp.test - flm_pred_test)^2)/mean(abs(resp.test)) # normalizing by the mean of the absolute value of y
    flm_mse_mean_sqr = mean((resp.test - flm_pred_test)^2)/mean(resp.test^2) # normalizing by the mean of the absolute value of y
    flm_mspe_lst = list("flm_mspe" = flm_mspe, "flm_mse_perc" = flm_mse_perc, "flm_mse_var" = flm_mse_var,
                        "flm_mean_abs" = flm_mse_mean_abs, "flm_mse_mean_sqr" = flm_mse_mean_sqr, "flm_run_time" = flm_run_time, 
                        "flm_train_MSE" = flm_train_mse)
    comp_mod_mspe[["flm"]] = flm_mspe_lst
  }

  ######## Poisson FLM ########
  if("gflm" %in% models){
    # Combining the training nonfd data with the FPC scores
    pois_pred_train = cbind(lin_var.train, std.scores_train)
    pois_flm_fit_start = unclass(Sys.time())
    # Poisson ridge regression which performs cross-validation for the penalty parameter.
    # Using alpha = 0 corresponds to the ridge penalty
    pois_flm = cv.glmnet(pois_pred_train, resp.train, family = "poisson", alpha = 0)
    t_pois_flm = unclass(Sys.time()) - pois_flm_fit_start
    pois_flm_run_time = t_pois_flm + t_fpca
    # Calculating training MSE
    pois_pred_train = predict(pois_flm, newx = pois_pred_train, type = "response", s = "lambda.min")
    pois_flm_train_mse = mean((resp.train-pois_pred_train)^2)
    # The test data
    pois_pred_test = cbind(lin_var.test, std.scores_test)
    # Prediction using the optimal lambda that minimizes MSPE
    pred_pois_flm = predict(pois_flm, newx = pois_pred_test, type = "response", s = "lambda.min")
    # Calculating MSPE's
    pois_flm_mspe = mean((resp.test-pred_pois_flm)^2)
    pois_flm_mse_perc = mean((resp.test - pred_pois_flm)^2/(resp.test)^2) # normalizing by the square y
    pois_flm_mse_var = mean((resp.test - pred_pois_flm)^2)/var(resp.test) # normalizing by the variance of y
    pois_flm_mse_mean_abs = mean((resp.test - pred_pois_flm)^2)/mean(abs(resp.test)) # normalizing by the mean of the absolute value of y
    pois_flm_mse_mean_sqr = mean((resp.test - pred_pois_flm)^2)/mean(resp.test^2) # normalizing by the mean of the absolute value of y
    pois_flm_mspe_lst = list("pois_flm_mspe" = pois_flm_mspe, "pois_flm_mse_perc" = pois_flm_mse_perc, "pois_flm_mse_var" = pois_flm_mse_var,
                        "pois_flm_mean_abs" = pois_flm_mse_mean_abs, "pois_flm_mse_mean_sqr" = pois_flm_mse_mean_sqr, 
                        "pois_flm_train_MSE" = pois_flm_train_mse,"pois_flm_run_time" = pois_flm_run_time)
    comp_mod_mspe[["gflm"]] = pois_flm_mspe_lst
  }
  return(comp_mod_mspe)
}





















