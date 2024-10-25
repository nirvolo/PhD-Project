model_comp_qr = function(train.lst_fd, test.lst_fd, nonfd.train, nonfd.test, fpca_thresh = 0.95, models, num_harm = 5,
                         train_quant = NULL, test_quant = NULL, sim, tau_val, cntr_fxns = T, pred_vars, DE_MEAN_RESP = F, algo = "br"){
 n_cnty = length(train.lst_fd) # number of counties in train and test
 n_knots = nrow(train.lst_fd[[1]]$coefs[,,1])
 # Initializing the dataframes by adding the first counties coefficients
 x1_train = train.lst_fd[[1]]$coefs[,,1]
 x2_train = train.lst_fd[[1]]$coefs[,,2]
 x1_test = test.lst_fd[[1]]$coefs[,,1]
 x2_test = test.lst_fd[[1]]$coefs[,,2]
 
 ##### Organize the excluded scores and organizing x1 and x2 into matrices ######s
 if(!sim){
   exld_scores = train.lst_fd[[1]]$excld_yrs_frm_trn
   shift_by = ncol(x1_train)
 }
 for(i in 2:n_cnty){
   fd_cnty_train = train.lst_fd[[i]]$coefs
   fd_cnty_test = test.lst_fd[[i]]$coefs
   x1_train = cbind(x1_train, fd_cnty_train[,,1])
   x2_train = cbind(x2_train, fd_cnty_train[,,2])
   x1_test = cbind(x1_test, fd_cnty_test[,,1])
   x2_test = cbind(x2_test, fd_cnty_test[,,2])
   # Adding the excluded scores for each county into the vector exld scores
   if(!sim){
     exld_scores = append(exld_scores, (train.lst_fd[[i]]$excld_yrs_frm_trn + shift_by))
     shift_by = ncol(x1_train)
   }
 }
 
 # Adding the x1 and x2 matrices to the train and test arrays
 # Arrays for to store the coefficient matrices for train and test
 arry_train = array(dim = c(n_knots, ncol(x1_train), 2))
 arry_test = array(dim = c(n_knots, ncol(x1_test), 2))
 arry_train[,,1] = x1_train
 arry_train[,,2] = x2_train
 arry_test[,,1] = x1_test
 arry_test[,,2] = x2_test
 # Creating the combined train and test fd object for fpca
 basis_fxns = train.lst_fd[[1]]$basis
 fd.train = fd(arry_train, basis_fxns)
 fd.test = fd(arry_test, basis_fxns)
 
 #### Perform FPCA with training data
 fpca_train = joint.fpca(fd.train, nharm = num_harm, thresh = fpca_thresh, centerfns = cntr_fxns, simulation = sim)
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
 if(!is.null(num_harm)){
   num_used_harm = num_harm
 } else if(is.null(ncol(fpc.scores_train))){
   #print("num scores is 1") 
   num_used_harm = 1
 } else{
   #print("I've entered ncol")
   num_used_harm = ncol(fpc.scores_train)
 }
 
 ####### Organizing the data for model fitting
 if(!sim){
   # Organizing the scalar data and the pc scores into one matrix
   # Adding an if statement to check the number of harmonics. If it's zero, don't include the fpc scores in train_pred
   if(!is.null(num_harm)){ 
     if(num_harm == 0) train_pred = as.matrix(cbind(1,nonfd.train[,pred_vars]))
     # The else to account for when nharm is greater than 0 but not NULL
     else train_pred = as.matrix(cbind(1,nonfd.train[,pred_vars], fpc.scores_train))
   } else{
     train_pred = as.matrix(cbind(1, nonfd.train[,pred_vars], fpc.scores_train)) # the data for the model
   }
   # the response variable
   if (!DE_MEAN_RESP) train_y = nonfd.train[,c("Yield")]
   else train_y = nonfd.train[,c("de_meaned_Yield")]
 } else{ # for simulated data
   train_pred = as.matrix(cbind(1,nonfd.train[,c("z")], fpc.scores_train))
   # the response variable
   train_y = nonfd.train[,c("y")]
 }
 
 # Combining the pred data with the response to give to the rq function
 train_data = data.frame(cbind(train_y, train_pred)) # rq function requires a data frame object
 colnames(train_data) = c("y", seq(1,ncol(train_pred)))
 
 ####### Fitting the functional quantile regression model
 train_qr_fit = rq(y ~ 0 + ., data = train_data, tau = tau_val, method = algo)
 
 ####### Calculate the estimated fpc scores for the test data and constructing the test data matrix for prediction #########
 pc_scores.test = proj.mfpca_nv(fpca_train$fpca, fd.test, cntr_fxns)$pc_scores[,1:num_used_harm]
 # Adding an if condition for simulation since the variable names are different
 if(!sim){
   # Data matrix for spatial VCM
   # Adding an if statement to check if the number of harmonics is zero. If it's zero, don't include scores
   if(num_used_harm == 0){
     test_X = as.matrix(cbind(1, nonfd.test[,pred_vars]))
   } else{
     test_X = as.matrix(cbind(1, nonfd.test[,pred_vars], pc_scores.test))
   }
   # Add a de_meaned_Yield test_y_demeaned
   if (!DE_MEAN_RESP) test_y = nonfd.test[,c("Yield")]
   else test_y = nonfd.test[,c("de_meaned_Yield")]
 } else{ # for simulated data
   test_X = as.matrix(cbind(1, nonfd.test[,c("z")], pc_scores.test))
   # the response variable
   test_y = nonfd.test[,c("y")]
 }

 # Changing the column names of train_X and test_X so they can be used for prediction.
 # This is because the coef names for the fit object must match the names of the dataset.
 # Also changing to a dataframe object since this is required for rq.
 train_pred = data.frame(train_pred)
 colnames(train_pred) = colnames(train_data)[-1]
 test_X = data.frame(test_X)
 colnames(test_X) = colnames(train_data)[-1]
 # Predicting train and test quantiles
 ypred_train = predict(train_qr_fit, train_pred) 
 ypred_test = predict(train_qr_fit, test_X)
 
 # Calculating the train and test MSE 
 if(!sim){ # for real data
   train_mse = sep(true_y = train_y, pred_y = ypred_train, tau = tau_val)
   test_mse = sep(true_y = test_y, pred_y = ypred_test, tau = tau_val)
 } else{
   train_mse = mean((train_quant-ypred_train)^2)
   test_mse = mean((test_quant-ypred_test)^2)
 }
 
 return(list("train_mse" = train_mse, "test_mse" = test_mse, "num_used_harm" = num_used_harm, "y_diff_train" = train_y-ypred_train, 
             "y_diff_test" = test_y-ypred_test, "train_fit" = train_qr_fit))
}

















