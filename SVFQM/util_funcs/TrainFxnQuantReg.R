# install.packages("quantreg")
# install.packages("rqPen")
# install.packages("hqreg")
# install.packges("expm")
# install.packages("data.table")
library(quantreg)
library(rqPen)
library(hqreg) # the package hqreg contains the function "hqreg_raw"
# this package is for the sqrtm function. If I don't run this, uses another sqrtm from a different package
# and it gives a numeric issue and recipricol condition number
library(expm) 
library(data.table) # this package is used by rqPen but the rqPen package doesn't load it for some reason

train_fxn_qr = function(train_fd, train_nonfd, full_fd, use_full_fd4mean = F, sim = F,
                           nharm = NULL, thresh = 0.80, d = 3, r = 1, 
                           spat_tri, cntr_fxns = F, DE_MEAN_RESP = T, pred_vars,
                           true_quant = NULL, tau, pnlty = "Ridge", n_lam = 10, train = T, scl = T, tmpmax_pars,lam_seq = NULL){
  if (length(tau) > 1) stop ("Error - for this function you must use a single tau value")
  # This function takes train_fd objects that contains county/year for which there is no data in nonfd
  # it uses all of it to find principal components, but take the scores and do training only on those that exist
  # in nonfd. The observations for county/year that are not in nonfd data would have $county NULL
  # sp_locs: All locations in the spatial domain
  # train: this variable is for predicting the quantiles which is necessary since calculation is different 
  # for train and test
  
  ###### Calculate the location-specific means of the training data #######
  fd_cntyIds = c()
  for (i in 1:length(train_fd)){  # iterate over counties in train fd
    fd_cntyIds = append(fd_cntyIds, train_fd[[i]]$CountyI)
  }
  
  train_fd_countiesId = unique(fd_cntyIds)
  
  if(!sim){
    train_counties = unique(train_nonfd$county)
  }
  
  n_knots = nrow(train_fd[[1]]$coefs[,,1]) # the number of knots for the smoothed data
  # Creating matrices to store the county means for the min and max temperatures
  # so that they can be used for the testing data later on.
  min.mean_mat = data.frame(matrix(0, nrow = n_knots, ncol = length(train_fd_countiesId)))
  max.mean_mat = data.frame(matrix(0, nrow = n_knots, ncol = length(train_fd_countiesId)))
  # Setting the names of the columns to the counties
  colnames(min.mean_mat) = train_fd_countiesId
  colnames(max.mean_mat) = train_fd_countiesId
  
  if(!sim){
    # Centering each county's coefficients by their means
    for(i in 1:length(train_fd)){
      train_cnty_data = train_fd[[i]] # the ith county data
      train_cntyId = train_cnty_data$CountyI
      # calculating the row means for the temperature which is the mean for each
      # of the spline coefficients
      cnty_min.mean = apply(train_cnty_data$coefs[,,1], 1, mean)
      cnty_max.mean = apply(train_cnty_data$coefs[,,2], 1, mean)
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
  # Creating the fd object that will be used in the FPCA
  fd_basis_fxns = train_fd[[1]]$basis # the basis fxns are the same for all counties
  joint_fd_obj = fd(comb_arry, fd_basis_fxns) # creates the "fd" object using the array of coefficients and basis functions.
  ##### Perform FPCA with the training data ####
  # Joint FPCA of the smoothed functional data X1 and X2. Since the data is already centered
  # above, I need to tell the joint.fpca function not to center the data which is why centerfns = F.
  fpca_train = joint.fpca(joint_fd_obj, nharm = nharm, thresh = thresh, centerfns = cntr_fxns, simulation = sim)

  if(!sim){
    # Including an if statement to check if there are excluding scores. 
    # If the length(exld_scores) == 0, then the command below returns an empty matrix.
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
    num_used_harm = 1
  } else{
    num_used_harm = ncol(fpc.scores_train)
  }
  
  ######## Calculating the penalty matrix Omega for the model fitting using its square root decomposition ###########
  train_basis_full = basis(spat_tri$V, spat_tri$Tr, d = d, r = r, as.matrix(train_nonfd[,c("long","lat")]))
  Q2_mat = train_basis_full$Q2
  B_mat = train_basis_full$B
  # Adding a check for dimensions 
  if(ncol(B_mat) != nrow(Q2_mat)) print("ERROR: dim of B is not the same as dim of Q2")
  Bstar_mat = B_mat%*%Q2_mat
  K_mat = train_basis_full$K
  if(nrow(Q2_mat) != nrow(K_mat)) print("ERROR: dim of B is not the same as dim of Q2")
  omega_pnlty = t(Q2_mat)%*%K_mat%*%Q2_mat
  # To make sure  we use the correct sqrtm function
  my_sqrtm <- expm::sqrtm
  # Changing type to numeric since sqrtm returns complex numbers but in my case it's equivalent to integers
  omega_sqrt = matrix(as.numeric(my_sqrtm(omega_pnlty)), nr = nrow(omega_pnlty), nc = ncol(omega_pnlty))
  if(!(max(abs(omega_sqrt%*%omega_sqrt - omega_pnlty)) < 1e-7)) print("ERROR: Square root of Omega is wrong")
  # Checking if the transpose of omega_sqrt is equal to itself
  if(!(max(abs(omega_sqrt - t(omega_sqrt))) < 1e-7)) print("ERROR: Omega not equal Omega^T")
  omega_sqrt_inv = solve(omega_sqrt, tol = 1e-21)
  if(!(max(abs(omega_sqrt_inv%*%omega_sqrt - diag(1, nrow = nrow(omega_sqrt), ncol = ncol(omega_sqrt)))) < 1e-7)){
    print(max(abs(omega_sqrt_inv%*%omega_sqrt - diag(1, nrow = nrow(omega_sqrt), ncol = ncol(omega_sqrt)))))
    print("ERROR: Inverse of Square Root of Omega is incorrect")
  }
  #### Organizing the data ####
  # Adding an if condition for simulation since the variable names are different
  if(!sim){
    # Organizing the scalar data and the pc scores into one matrix
    # Adding an if statement to check the number of harmonics. If it's zero, don't include the fpc scores in train_pred
    if(!is.null(nharm)){ 
      if(nharm == 0) train_pred = as.matrix(cbind(1,train_nonfd[,pred_vars]))
      # The else to account for when nharm is greater than 0 but not NULL
      else train_pred = as.matrix(cbind(1,train_nonfd[,pred_vars], fpc.scores_train))
    } else{
      train_pred = as.matrix(cbind(1, train_nonfd[,pred_vars], fpc.scores_train)) # the data for the model
    }
    # the response variable
    if (!DE_MEAN_RESP) train_resp = train_nonfd[,c("Yield")]
    else train_resp = train_nonfd[,c("de_meaned_Yield")]
    # Calculating the reparametrized data matrix using the penalty matrix Omega.
    # For the rq.pen.cv function, the data matrix type is required to be a numeric matrix.
    X_tilda = as.matrix(kr(train_pred,as.matrix(Bstar_mat)%*%omega_sqrt_inv, byrow = T))
  } else{
    # for simulated data
    train_pred = as.matrix(cbind(1,train_nonfd[,c("z")], fpc.scores_train))
    # the response variable
    train_resp = train_nonfd[,c("y")]
    # Calculating the reparametrized data matrix using the penalty matrix Omega
    X_tilda = as.matrix(kr(train_pred,as.matrix(Bstar_mat)%*%omega_sqrt_inv, byrow = T))
  }
  
  ###### Fit the Quantile Regression Model using the predictor variables above ######
  # intrcpt is set to FALSE here so that the transformation of the coefficients is done correctly
  train_fit = rq.pen.cv_new(X_tilda, train_resp, penalty = pnlty, nlambda = n_lam, tau = tau, scalex = scl, intrcpt = F, lambda = lam_seq)
  # Obtaining the lambda index corresponding to the lambda that minimimizes the CV error
  min_cv_lam_ind = train_fit$btr$lambdaIndex
  # Extracting the training coefficients. The first index corresponds to acccessing the list
  # and the second is the coefs. The reason for this is because the model fit is labeled
  # by "tau0.5a" where 0.5 is the quantile and this will change so not general enough for other values.
  train_coefs = train_fit$fit$models[[1]][[1]][,min_cv_lam_ind]
  ##############################################
  # Calculate the coefficient matrix
  nr_omega = nrow(omega_sqrt_inv)
  n_coefs = length(train_coefs)
  theta_tilda = matrix(train_coefs, nr = nr_omega, nc = n_coefs/nr_omega)
  if(ncol(omega_sqrt_inv) != nrow(theta_tilda)) print("ERROR") 
  theta_mat = omega_sqrt_inv%*%theta_tilda
  # Calculating the predicted quantiles
  pred_otpt = predict_svqfm(spline_coefs = theta_mat, Xpred = train_pred,
                            triang = spat_tri, B = B_mat, Q2 = Q2_mat, train = train)
  train_yhat = pred_otpt$yhat
  train_est_etas = pred_otpt$eta
  # Calculating the train MSE which compares the estimated quantiles y_hat to the true.
  # This is only the case for simulations, otherwise calculate the MSE using the check loss
  # function
  if(!sim){
    # Calcuating the training MSE
    train_mse = sep(true_y = train_resp, pred_y = train_yhat, tau = tau)
  } else{
    # Calculating the training MSE
    train_mse = mean((train_yhat - true_quant)^2)
  }
  # Including the response and the train_pred in the output list so this can be used for the comparison model
  train_otpt = list("fpca_obj" = fpca_train$fpca, "cnty_means" = cnty_means, 
                    "train_mse" = train_mse, "num_used_harm" = num_used_harm, "spline_coefs" = theta_mat,
                    "train_B" = B_mat, "train_Q2" = Q2_mat, "train_y" = train_resp, "train_pred" = train_pred, "train_etas" = train_est_etas,
                    "y_diff_train" = train_resp-train_yhat, "train_fit" = train_fit, "train_yhat" = train_yhat)
  return(train_otpt)
}





