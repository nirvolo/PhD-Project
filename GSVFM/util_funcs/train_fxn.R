train_fxn_allfd = function(train_fd, train_nonfd, full_fd, use_full_fd4mean = F, sim = F,
                           nharm = NULL, thresh = 0.80, d = 3, r = 1, lambda = exp(seq(log(0.00001), log(100), length.out=10)), 
                           family = gaussian(), sp_tri, sp_locs = predictor_vars[,c("long","lat")], B0_mat = alldata_B0_mat, cntr_fxns = F, 
                           coef_comp = "complex", DE_MEAN_RESP = T, pred_vars = pred_vars){
  # This function takes train_fd objects that contains county/year for which there is no data in nonfd
  # it uses all of it to find principal components, but take the scores and do training only on those that exist
  # in nonfd. The observations for county/year that are not in nonfd data would have $county NULL
  # sp_locs: All locations in the spatial domain
  
  ###### Calculate the location-specific means of the training data #######
  fd_cntyIds = c()
  for (i in 1:length(train_fd)){  # iterate over counties in train fd
    fd_cntyIds = append(fd_cntyIds, train_fd[[i]]$CountyI)
  }
  # The unique county IDs
  train_fd_countiesId = unique(fd_cntyIds)
  n_knots = nrow(train_fd[[1]]$coefs[,,1]) # the number of knots for the smoothed data
  # Creating matrices to store the county means for the min and max temperatures
  # so that they can be used for the testing data later on.
  min.mean_mat = data.frame(matrix(0, nrow = n_knots, ncol = length(train_fd_countiesId)))
  max.mean_mat = data.frame(matrix(0, nrow = n_knots, ncol = length(train_fd_countiesId)))
  # Setting the names of the columns to the counties
  colnames(min.mean_mat) = train_fd_countiesId
  colnames(max.mean_mat) = train_fd_countiesId
  
  t_start = unclass(Sys.time()) # Unclass to get time as double in seconds
  if(!sim){
    # Centering each county's coefficients by their means
    for(i in 1:length(train_fd)){
      train_cnty_data = train_fd[[i]] # the ith county data
      train_cntyId = train_cnty_data$CountyI
      # calculating the row means for the temperature which is the mean for each
      # of the spline coefficients
      if(use_full_fd4mean){
        cnty_min.mean = apply(full_fd[[i]]$coefs[,,1], 1, mean)
        cnty_max.mean = apply(full_fd[[i]]$coefs[,,2], 1, mean)
      } else{  # Mean per county
        # There are some counties with only 1 year of data which causes an issue
        # with the apply function. These kind of counties would have county name NA
        # because they are not in regdat but we still use them for FPCA. They don't 
        # make it through the filtering in ag_cv fxn since we use the num_obs from regdat.
        if(is.null(dim(train_cnty_data$coefs[,,1]))){
          cnty_min.mean = train_cnty_data$coefs[,,1]
          cnty_max.mean = train_cnty_data$coefs[,,2]
        } else{
          cnty_min.mean = apply(train_cnty_data$coefs[,,1], 1, mean)
          cnty_max.mean = apply(train_cnty_data$coefs[,,2], 1, mean)
        }
      }
      # Centering the spline coefficients by their location-specific mean
      train_fd[[i]]$coefs[,,1] = train_cnty_data$coefs[,,1] - cnty_min.mean
      train_fd[[i]]$coefs[,,2] = train_cnty_data$coefs[,,2] - cnty_max.mean
      # Adding the means to the matrix
      min.mean_mat[,as.character(train_cntyId)] = cnty_min.mean
      max.mean_mat[,as.character(train_cntyId)] = cnty_max.mean
    }
  }

  # Creating a list to store the min and max mean matrices
  cnty_means = list()
  cnty_means[[1]] = min.mean_mat
  cnty_means[[2]] = max.mean_mat
  
  # Creating matrices to store the min and max coefficients.
  # Initializing the matrix by adding the first counties coefficients
  min.temp_mat = train_fd[[1]]$coefs[,,1]  # after demean
  max.temp_mat = train_fd[[1]]$coefs[,,2]
  # Each column of coefs matrix correspond to a year in this county
  # Use the list of years to exclude to later omit these columns from the scores
  # given to fit.gsvcm. For each element of the list, take the exld indices
  # and keep the number of columns in min_temp_mat so far and use to shift the
  # indices of the next list to exlude
  if(!sim){
    exld_scores = train_fd[[1]]$excld_yrs_frm_trn
    shift_by = ncol(min.temp_mat)
  }
  for(j in 2:length(train_fd_countiesId)){
    min.temp_mat = cbind(min.temp_mat, train_fd[[j]]$coefs[,,1])
    max.temp_mat = cbind(max.temp_mat, train_fd[[j]]$coefs[,,2])
    if(!sim){
      exld_scores = append(exld_scores, (train_fd[[j]]$excld_yrs_frm_trn + shift_by))
      shift_by = ncol(min.temp_mat)
    }
  }
  ##### Combining the data from ALL counties for FPCA #####
  comb_arry = array(dim = c(n_knots, dim(min.temp_mat)[2], 2))
  # Adding the matrices to the array
  comb_arry[,,1] = min.temp_mat
  comb_arry[,,2] = max.temp_mat
  t_cntring = unclass(Sys.time()) - t_start
  t_start_fpca = unclass(Sys.time())
  # Creating the fd object that will be used in the FPCA
  fd_basis_fxns = train_fd[[1]]$basis # the basis fxns are the same for all counties
  joint_fd_obj = fd(comb_arry, fd_basis_fxns) # creates the "fd" object using the array of coefficients and basis functions.
  
  ##### Perform FPCA with the training data ####
  # Joint FPCA of the smoothed functional data X1 and X2. Since the data is already centered
  # above, I need to tell the joint.fpca function not to center the data which is why centerfns = F.
  fpca_train = joint.fpca(joint_fd_obj, nharm = nharm, thresh = thresh, centerfns = cntr_fxns, simulation = sim, scl = scl_scores)
  if(!sim){
    # Including an if statement to check if there are excluding scores. 
    if(length(exld_scores) > 0){
        # The if statement is to account for the case where there is only 1 FPC score so the scores is a vector not a matrix
        if(is.null(dim(fpca_train$pc_scores))) fpc.scores_train = fpca_train$pc_scores[-exld_scores]
        else fpc.scores_train = fpca_train$pc_scores[-exld_scores,]  # Exclude scores for obs with no nonfd data
    } else{
        fpc.scores_train = fpca_train$pc_scores # Exclude scores for obs with no nonfd data
    }
  } else{
      fpc.scores_train = fpca_train$pc_scores # Exclude scores for obs with no nonfd data
  }
  # The first condition is only for the case where I give the function a set
  # number of harmonics, 0,1,2,3..... Otherwise it'll select the nharm based
  # on the number of columns in fpc.scores_train. The default for nharm is null,
  # only gets an argument when I don't want to use threshold.
  if(!is.null(nharm)){
    num_used_harm = nharm
  } else if(is.null(ncol(fpc.scores_train))){
    print("num scores is 1") 
    num_used_harm = 1
  } else{
    num_used_harm = ncol(fpc.scores_train)
  }
  t_fpca = unclass(Sys.time()) - t_start_fpca
  
  #### Organizing the data ####
  # Adding an if condition for simulation since the variable names are different
  if(!sim){
    # Organizing the scalar data and the pc scores into one matrix
    # Adding an if statement to check the number of harmonics. If it's zero, don't include the fpc scores in train_pred
     if(!is.null(nharm)){ 
      if(nharm == 0) train_pred = as.matrix(train_nonfd[,pred_vars])
      # The else to account for when nharm is greater than 0 but not NULL
      else train_pred = as.matrix(cbind(train_nonfd[,pred_vars], fpc.scores_train))
    } else{
      train_pred = as.matrix(cbind(train_nonfd[,pred_vars], fpc.scores_train)) # the data for the model
    }
    # the response variable
    if (!DE_MEAN_RESP) train_resp = train_nonfd[,c("Yield")]
    else train_resp = train_nonfd[,c("de_meaned_Yield")]
  } else{
    # for simulated data
    train_pred = as.matrix(cbind(train_nonfd[,c("z")], fpc.scores_train))
    # the response variable
    train_resp = train_nonfd[,c("y")]
  }
  # Adding the intercept
  train_pred = as.matrix(cbind(1, train_pred))
  # The coordinates of the training locations
  train_locs = as.matrix(train_nonfd[,c("long","lat")])
  
  ###### Fit the SVCM using the predictor variables above ######
  # Measuring the run time for the training model
  t_start_train = unclass(Sys.time())
  train_fit = fit.gsvcm(train_resp, train_pred, train_locs, sp_tri$V, sp_tri$Tr,
                        d, r, lambda, family, off = 0)
  t_train = unclass(Sys.time()) - t_start_train
  train_run_time = t_cntring + t_fpca + t_train
  # Calculating MSE for the training set
  train_yhat = predict(train_fit, train_pred, as.matrix(train_nonfd[,c("long","lat")]))
  train_mse = mean((train_yhat - train_resp)^2)
  
  # List of outputs for the function
  train_otpt_list = list("fpca_obj" = fpca_train$fpca, "fpc_scores" = fpc.scores_train, "fd_obj" = joint_fd_obj,
                         "svfm_fit" = train_fit, "cnty_means" = cnty_means, "train_mse" = train_mse, "num_used_harm" = num_used_harm,
                         "svfm_run_time" = train_run_time, "fpca_eigenfuncs" = fpca_train$fpca, "centered_train_fd" = train_fd)
  return(train_otpt_list)
}
















