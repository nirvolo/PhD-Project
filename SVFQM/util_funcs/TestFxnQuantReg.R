test_fxn_qr = function(test_fd, test_nonfd, cnty_means, train_fpca_obj, train_beta_mat, num_harm, cntr_fxns = F, sim = F, DE_MEAN_RESP = T, 
                       pred_vars = pred_vars, true_quant = NULL, train = F, spline_coefs, train_B, train_Q2, triang, tau){
  ####### Center the testing data by their location specific means ########
  # The means for min and max, respectively
  min_means = cnty_means[[1]]
  max_means = cnty_means[[2]]
  if(!sim){
    # Iterating over all the counties in test_fd
    for(i in 1:length(test_fd)){
      # The ith county
      test_cntyId = test_fd[[i]]$CountyI
      # Centering the data
      test_fd[[i]]$coefs[,,1] = test_fd[[i]]$coefs[,,1] - unlist(min_means[, as.character(test_cntyId)]) # for min temp - unlist convert to numeric vector (double)
      test_fd[[i]]$coefs[,,2] = test_fd[[i]]$coefs[,,2] - unlist(max_means[, as.character(test_cntyId)]) # for max temp
    }
  }
  # Combining the data from all the testing counties
  # Checking if the dimension of the test data is NULL since there might be one observation
  # for cross-validation and you can't take the dimension of a vector
  if(is.null(dim(test_fd[[1]]$coefs[,,1]))){
    n_knots = length(test_fd[[1]]$coefs[,,1]) # the number of knots for the smoothed data
  } else{
    n_knots = nrow(test_fd[[1]]$coefs[,,1]) # the number of knots for the smoothed data
  }
  comb_test_arry = array(dim = c(n_knots, dim(test_nonfd)[1], 2))
  # Creating matrices to store the min and max coefficients.
  # Initializing the matrix by adding the first counties coefficients
  min.temp_test = test_fd[[1]]$coefs[,,1]
  max.temp_test =  test_fd[[1]]$coefs[,,2]
  for(j in 2:length(test_fd)){
    # Combining the coefficients for the min and max temp.
    min.temp_test = cbind(min.temp_test, test_fd[[j]]$coefs[,,1])
    max.temp_test = cbind(max.temp_test, test_fd[[j]]$coefs[,,2])
  }
  # Adding the matrices to the array
  comb_test_arry[,,1] = min.temp_test
  comb_test_arry[,,2] = max.temp_test
  
  # Creating the fd object that will be used in the FPCA
  fd.test_basis_fxns = test_fd[[1]]$basis # the basis fxns are the same for all counties
  joint.test_fd_obj = fd(comb_test_arry, fd.test_basis_fxns) # creates the "fd" object using the array of coefficients and basis functions.
  
  ####### Calculate the estimated fpc scores for the test data #########
  pc_scores.test = proj.mfpca_nv(train_fpca_obj, joint.test_fd_obj, cntr_fxns)$pc_scores[,1:num_harm]

  ####### Constructing the test data matrix for prediction ########
  # Adding an if condition for simulation since the variable names are different
  if(!sim){
    # Data matrix for spatial VCM
    # Adding an if statement to check if the number of harmonics is zero. If it's zero, don't include scores
    if(num_harm == 0){
      test_X = as.matrix(cbind(1, test_nonfd[,pred_vars]))
    } else{
      test_X = as.matrix(cbind(1, test_nonfd[,pred_vars], pc_scores.test))
    }
    # Add a de_meaned_Yield test_y_demeaned
    if (!DE_MEAN_RESP) test_y = test_nonfd[,c("Yield")]
    else test_y = test_nonfd[,c("de_meaned_Yield")]
  } else{
    # for simulated data
    test_X = as.matrix(cbind(1, test_nonfd[,c("z")], pc_scores.test))
    # the response variable
    test_y = test_nonfd[,c("y")]
  }
  
  ####### Predict the response variable or all values of the test data and calculate MSE ######
  # Predict the estimated test quantiles
  test_otpt = predict_svqfm(spline_coefs = spline_coefs, Xpred = test_X, Spred = test_nonfd[,c("long","lat")], 
                            triang = triang, Q2 = train_Q2, train = train)
  test_yhat = test_otpt$yhat
  
  # MSE for the predicted y for test data. MSE is different depending on whether it's a simulation or not
  if(!sim){
    # MSE for real data
    test_mse = sep(true_y = test_y, pred_y = test_yhat, tau = tau)
  } else{
    # MSE for simulation data
    num_pred = length(test_yhat)
    sum_y_great_est = sum(test_y > test_yhat)
    test_mse = (sum_y_great_est/num_pred) - (1-tau)
  }
  # Including the response and the train_pred in the output list so this can be used for the comparison model
  return(list("test_MSE" = test_mse, "test_y" = test_y, "test_pred" = test_X))
}

