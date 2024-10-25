ag.data_cv.fxn_nopresmooth = function(fundat_all, nonfd, n_fold, n_yr = 13, thresh = 0.8, sp_tri = kns_tri2, iter = 1, deg = 3, 
                                      DE_MEAN_RESP = T, pred_vars, reval, log.lam = c(-15.5,5), fld_lim = 5, nharm = NULL){
  # fundat_all: raw functional data after preprocessing but without smoothing, not an fd objects list
  #             Contains functional data for county/years that are not in nonfd - these are at the end of fundat_all
  set.seed(12876+iter*123)
  cat("starting", n_fold, "folds iteration", iter, "\n")
  cat("=========================================", "\n")
  # Percentage of data to use for training
  ratio = 1 - (1/n_fold)
  
  # Vector to store the train MSE, test MSPE and train run time for each fold
  mspe = c()
  train_runtime = c()
  svfm_train_mse = c() # I only calculate the MSE for training
  svfm_mspe = data.frame(matrix(NA, nr = n_fold, nc = 5))
  colnames(svfm_mspe) = c("test_MSE", "test_mse_perc", "test_mse_var",
                          "test_mse_mean_abs", "test_mse_mean_sqr")
  for(fld in 1:fld_lim){
  #for(fld in 1:1){
    cat("Fold", fld, "\n")
    cat('------------', "\n")
    tst_ids = c()
    # Vector to store cnty ids that have less than a certain number of years
    excl_cnty_ids = c()
    
    for(cnty_id in unique(nonfd$CountyI)){  
      n_obs = nrow(nonfd[nonfd$CountyI == cnty_id,])
      # Not including counties that have less than 5 years of data
      if(n_obs < 5){ 
        excl_cnty_ids = append(excl_cnty_ids, cnty_id)
        next
      }
      # Number of training and testing observations for this county
      n_train = floor(ratio*n_obs)
      n_test = n_obs - n_train
      # The train and test indices for the current county
      # In each fold we will get a different set of indices for training with slight
      # probability for repetition.
      trn_inds = sample(1:n_obs, n_train)
      tst_inds = (1:n_obs)[-trn_inds]
      tst_ids = append(tst_ids, nonfd[nonfd$CountyI == cnty_id,][tst_inds,]$Id)  # Add obs Id to tst_ids
    } # end of for loop for counties  
    
    ########## Organizing the non-functional data ##############
    # The nonfd data for the current county and window
    nonfd_test = nonfd[nonfd$Id %in% tst_ids,]
    tst_fd = fundat_all[fundat_all$Id %in% tst_ids,]
    # Excluding the counties that have less than 5 years of data in the training data 
    # so that there isn't a possibility of including counties with one year of data
    nonfd_train = nonfd[!nonfd$Id %in% tst_ids,]
    nonfd_train = nonfd_train[!nonfd_train$CountyI %in% excl_cnty_ids,]
    trn_fd = fundat_all[!fundat_all$Id %in% tst_ids,]
    trn_fd = trn_fd[!trn_fd$CountyI %in% excl_cnty_ids,]
    # Smoothing the functional data
    smoothed_trn_fd = fd_smooth_2steps_trn_tst(trn_fd, n_obs = dim(trn_fd)[1]/365, reval = reval, llam.lims = log.lam)
    min_opt_lam = smoothed_trn_fd$min_opt_lambda
    max_opt_lam = smoothed_trn_fd$max_opt_lambda
    comb_trn_fd_lst = smoothed_trn_fd$comb_fd
    # Using the optimal lambdas found in training set to smooth test set
    smoothed_tst_fd = fd_smooth_2steps_trn_tst(tst_fd, n_obs = dim(trn_fd)[1]/365, min_opt_lamb = min_opt_lam, max_opt_lamb = max_opt_lam, 
                                               reval = reval, llam.lims = log.lam)
    comb_tst_fd_lst = smoothed_tst_fd$comb_fd
    # Changing the column names of nonfd train and test
    colnames(nonfd_train) = colnames(nonfd)
    colnames(nonfd_test) = colnames(nonfd)
    
    ### Centering the yield data by the mean by year ###
    # Add column yld_mean_per_yr to train_nonfd
    nonfd_train$yld_mean_per_yr = ave(nonfd_train$Yield, nonfd_train$Year)
    # add De-Meaned column to train
    nonfd_train$de_meaned_Yield = nonfd_train$Yield - nonfd_train$yld_mean_per_yr
    # Add yld_mean_per_yr to nonfd_test
    nonfd_test$yld_mean_per_yr = NA
    for (yr in unique(nonfd_test$Year)){
      nonfd_test[nonfd_test$Year == yr,]$yld_mean_per_yr = unique(nonfd_train[nonfd_train$Year == yr,]$yld_mean_per_yr)
    }
    # Add column for De-meaned Yield in nonfd_test
    nonfd_test$de_meaned_Yield = nonfd_test$Yield - nonfd_test$yld_mean_per_yr
    
    # Fitting the training model
    #cat("The dimension of train_nonfd before training function is,",dim(na.omit(nonfd_train)),"\n")
    train_mod = train_fxn_allfd(train_fd = comb_trn_fd_lst, train_nonfd = na.omit(nonfd_train), full_fd = fundat_all, 
                                thresh = thresh, sp_tri = sp_tri, use_full_fd4mean = F, d = deg, DE_MEAN_RESP = DE_MEAN_RESP, pred_vars = pred_vars,
                                nharm = nharm)
    # Storing the training MSE for each fold
    svfm_train_mse = append(svfm_train_mse, train_mod$train_mse)
    cat("The train MSE is", train_mod$train_mse, "\n")
    mean_cntyIds = colnames(train_mod$cnty_means[[1]])
    nonfd_countyIds = unique(nonfd_train$CountyI)
    # Prediction using the testing data
    test_res = test_fxn_allfd(comb_tst_fd_lst, na.omit(nonfd_test), train_mod$cnty_means, train_mod$fpca_obj,
                                     train_mod$svfm_fit, num_harm = train_mod$num_used_harm, DE_MEAN_RESP = DE_MEAN_RESP, pred_vars = pred_vars)
    test_mspe_res = unlist(test_res)[1:5]
    svfm_mspe[fld,] = test_mspe_res
    mspe[fld] = svfm_mspe[fld,"test_MSE"]
    # Saving the train run time for the current fold
    train_runtime = train_mod$svfm_run_time
    cat("The MSPE for the current fold is", mspe[fld], "\n")
  } # End of for loop for folds
  # Storing the number of harmonics. I only need to do this once at the end of the folds since that same number of
  # harmonics is used for all folds.
  svfm_num_harm = train_mod$num_used_harm
  # Return all res so we can see normalized as well
  return(list("test_MSPE" = svfm_mspe, "train_runtime" = train_runtime, "train_MSE" = svfm_train_mse,"num_used_harm" = svfm_num_harm, "train_mod" = train_mod,
              "test_yhat" = test_res$test_yhat, "test_y" = test_res$test_y, "test_locs" = nonfd_test[,c("long","lat")]))
}







