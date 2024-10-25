# packages
# install.packages("spatstat.geom")
#install.packages("combinat")
#devtools::install_github('FIRST-Data-Lab/Triangulation') for installing triangulation package
library(spatstat.geom) # this is for the function cross_dist.
library(MASS) # for mvrnorm
library("fda")
library(combinat)
library('BPST') # This package is for the basis function
library('MGLM') 
library(plot3D)
library(devtools)
library(Triangulation)
library(tidyverse)
library(ggplot2)
library(raster)
library(RColorBrewer)
library(ggplot2)
library(classInt)
library(rnaturalearth)
library(sf)
library(rgdal)
library(rgeos)
library(tigris)
library(maps)
library(quantreg) # This is for the rq functions
#install.packages("hqreg")
library(quantreg)
library(rqPen)
library(hqreg) # the package hqreg contains the function "hqreg_raw"
library(expm) # this package is for the sqrtm function. If I don't run this, uses another sqrtm from a different package
              # and it gives a numeric issue and recipricol condition number
library(data.table) # this package is used by rqPen but the rqPen package doesn't load it for some reason


source("common_utils/preprocessing_function.R") 
source_rmd("common_utils/KansasMaps.Rmd")
source("SVFQM/util_funcs/TrainFxnQuantReg.R")
source("SVFQM/util_funcs/TestFxnQuantReg.R")
source("SVFQM/util_funcs/predict.svqfm.R")
source("SVFQM/util_funcs/ModelCompQR.R")
source("SVFQM/rqPen/rq.pen.cv.fxn.R")
source("SVFQM/rqPen/workHorse.R")

matern_corr = function(phi, d, kap = 1){
  if(d == 0) return(1)
  else return((1/(gamma(kap)*(2^(kap-1))))*((d/(phi/111))^kap)*besselK(d/(phi/111),kap))
}

eigen_basis = function(r, t){
  # Calculates the rth order eigenfunction for X1 and X2
  # r: the order of the eigenfunction
  # t: the time points
  max_t = t[length(t)]
  f1 = cos(2*pi*r*t/max_t)/sqrt(max_t)
  f2 = sin(2*pi*r*t/max_t)/sqrt(max_t)
  return(rbind(f1,f2))
}

coef_fxns = function(locs, coef){
  # Calculates the spatially varying coefficients given the spatial locations.
  # locs: spatial locations
  s1 = locs[,1] # longitude
  s2 = locs[,2] # latitude
  # The number of locations
  n_locs = nrow(locs)
  # coefficients alpha0,alpha1,beta1,...,beta5
  if(coef == "basic"){
    alpha0 = 1.5*s1 + 1.2*s2
    alpha1 = 2*s1 + 2*s2
    b1 = 4*s1 + 6*s2
    b2 = 3*s1 + 5*s2
    b3 = 3*s1 + 2*s1 
    b4 = 2*s1 + 5*s2
    b5 = 2*s1 + 2.5*s2
  } else if(coef == "complex"){
    alpha0 = 2*sin((pi*s1)/6)/4
    alpha1 = 2*sin(s1 + s2)/4
    b1 = (2/81)*(9-(3-s1)^2)*(9-(3-s2)^2)/11
    b2 = 2*sin((pi*(s1 + s2))/6)/4
    b3 = (1/2)*(8+(4-s1)^2)*(8-(4-s2)^2)/(660)
    b4 = sin(2*pi*(s1 + s2))/2
    b5 = (1/3)*(3+(2-sqrt(abs(s1)))^2)*(3-(2-sqrt(abs(s2)))^2)/8
  } else if(coef == "const"){
    alpha0 = rep(0.5, n_locs)
    alpha1 = rep(1,n_locs)
    b1 = rep(2,n_locs)
    b2 = rep(-2,n_locs)
    b3 = rep(-2.5,n_locs) 
    b4 = rep(-1.3,n_locs)
    b5 = rep(2.5,n_locs)
  }
  factor = 0.4/5 # constant case
  if(coef == "basic") factor = 0.1/6.5
  else if(coef == "complex") factor = 3/6
  
  return(matrix(factor*c(alpha0,alpha1,b1,b2,b3,b4,b5), nrow = dim(locs)[1], ncol = 7))
}

# Creating a function that simulates the true data and the data with noise
sim_fxn_qr = function(N, phi, v = sqrt(0.2), sigma = NULL, sp_locs, sp_bnd = kns_mat_scl, p = 5, years = 10, 
                       num_iter = 1, n_folds = 10, fld_lim = 10, n_pars = 7, coef_comp = "basic", nharm = 5, snr, 
                       comp_model = T, spat_corr, tri_fin, d = 3, r = 1, models, reval = c(0,1), log.lam = c(-15.5,-8), tau,
                       pnlty = "Ridge", n_lam = 10, qr_algo = "br", scl = T, tmpmax_pars = c(1e5,-12.1, 2.35)){
  # n: sample size
  # N: number of time points
  # phi: Matern covariance parameter
  # v: sd of error (functional predictor)
  # sigma: sd of error (response)
  # sp_locs: matrix of spatial locations
  # p: number of eigenfunctions and pc scores
  # years: number of observations per location
  # num_iter: the number of iterations
  # tau: the quantile we want to estimate
  # scl: whether or not to scale the data using rq.pen
  
  n_locs = dim(sp_locs)[1] # number of locations
  n = n_locs*years # number of observations
  # Creating a vector to store the County ID's. There are n_locs so there are n_locs County ID's
  cnty_ids = seq(1,n_locs)
  # Creating vectors to store the MSPE and run time for each models. These vectors will contain the values
  # for all replications and folds
  svqfm_train_mse = c()
  svqfm_mspe = c()
  fqr_train_mse = c()
  fqr_mspe = c()

  ######### Triangulation #########
  spat_tri = TriMesh(Pt = sp_bnd, n = tri_fin)
  
  ####################### Beginning of simulation for loop ##################
  for(iter in 1:num_iter){
    cat("Data generation", iter, "\n")
    set.seed(12876+iter*123)
    ########### Generate the Matern covariance matrix ###########
    # matrix of distances between all pairs of points
    dist_mat = crossdist(sp_locs[,1], sp_locs[,2], sp_locs[,1], sp_locs[,2])
    # Including an if statement here to check if the setting is to include
    # spatial covariance or no spatial covariance. 
    if(spat_corr){
      # Spatial covariance for all pairs of distaces
      spat_cov = apply(dist_mat, c(1,2), matern_corr, phi = phi)
    } else{
      # If there is no spatial covariance, the covariance matrix is 
      # simply an identity matrix.
      spat_cov = diag(nrow(dist_mat))
    }
    
    ######### Generating the eigenfunctions #################
    time_pts <- seq(reval[1], reval[2], len = N)
    # Generating the eigenfunctions by creating
    # a matrix where each row corresponds to the rth order
    # eigenfunction. pc1 is for X1 and pc2 is for X2.
    pc1 = matrix(0, nrow = p, ncol = N)
    pc2 = matrix(0, nrow = p, ncol = N)
    for(r in 1:p){
      pc_fxns = eigen_basis(r, t = time_pts)
      pc1[r,] = pc_fxns[1,]
      pc2[r,] = pc_fxns[2,]
    }
    
    ############### Generation of functional data X ##################### 
    # Mean of PC scores. The mean at each location is 0
    pc.score_mu = rep(0, n_locs)
    # Variances for the 5 principal component scores
    pc_var = c(9,4,3,2,1)
    # An empty list that will store the pc scores matrix
    # for the kth observation
    pc.scores_list = list()
    X01 <- data.frame(matrix(nr = n, nc = N)) # functional predictors with no measurement errors
    X02 <- data.frame(matrix(nr = n, nc = N)) # functional predictors with no measurement errors
    # Creating a dataframe for scoring the pc scores
    scores_df = data.frame(matrix(0, nr = n, nc = nharm))
    
    ######### Generating the principal component scores #########
    for(k in 1:years){ # Iterating over the number of years.
      pc_scores = matrix(0, nrow = p, ncol = n_locs) # matrix of all the pc scores
      # Simulating the p pc score vectors from MVN(0,var*Sigma)
      for(j in 1:p){
        pc_scores[j,] = mvrnorm(n = 1, mu = pc.score_mu, Sigma = pc_var[j]*spat_cov)
      }
      # Adding the kth pc scores matrix to the list
      pc.scores_list[[k]] = pc_scores
      ######### Generating the functional data X #############
      for(l in 1:n_locs){ # iterating over spatial locations
        scores = pc_scores[,l] # the pc scores for location l
        scores_df[(l-1)*years + k,] = scores
        for(t in 1:N){ # iterating over number of times
          # Calculating the elements X1_k(s_l;t) and X2_k(s_l;t)
          X01[(l-1)*years + k,t] =  scores%*%pc1[,t]
          X02[(l-1)*years + k,t] =  scores%*%pc2[,t]
        } 
        # Adding the county ID's to the functional data
        X01[(l-1)*years + k, "CountyI"] = cnty_ids[l]
        X02[(l-1)*years + k, "CountyI"] = cnty_ids[l]
        # Adding the long & lat to the functional data
        X01[(l-1)*years + k, c("long", "lat")] = sp_locs[l,]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        X02[(l-1)*years + k, c("long", "lat")] = sp_locs[l,]
      }
    }
    ## functional data with measurement errors N(0,0.2^2)
    # Need to add 3 columns to for county, long and lat
    X1 = data.frame(matrix(0, nrow = n, ncol = N)) 
    X2 = data.frame(matrix(0, nrow = n, ncol = N))
    fd_msr_error = data.frame(matrix(rnorm(n*N, sd=v), nr = n, nc = N))
    # Adding the measurement error and the county and location columns to X1 and X2
    X1[,1:N] <- X01[,1:N] + fd_msr_error
    X2[,1:N] <- X02[,1:N] + fd_msr_error
    X1[,c("CountyI", "long", "lat")] = X01[,c("CountyI", "long", "lat")]
    X2[,c("CountyI", "long", "lat")] = X02[,c("CountyI", "long", "lat")]
    # Adding a column for the Id component of the data
    X1[,"Id"] = seq(1, n)
    X2[,"Id"] = seq(1, n)
    # This code is changing the format of the X matrices from each column being a day to each row being a day.
    X1_new = X1 %>% pivot_longer(cols = colnames(X1)[1:N], names_to = "day", values_to = "x1")
    X2_new = X2 %>% pivot_longer(cols = colnames(X2)[1:N], names_to = "day", values_to = "x2")
    X_new = data.frame(cbind(X1_new, X2_new[,"x2"]))

    ############### Generating scalar predictors, coefficient functions, and response variable #######################
    # Creating the dataframe to combine all the non fd data
    sim_nonfd = data.frame(matrix(0, nr = n, nc = 6))
    colnames(sim_nonfd) = c("y", "z", "CountyI", "Id", "long", "lat")
    # Adding the Id variable to sim_nonfd
    sim_nonfd$Id = seq(1,n)
    # Generating the scalar predictor variable Z ~ U(0,3) where
    sim_nonfd$z = runif(n, min = 0, max = 3)
    
    ###### Generating the coefficients for the model ######
    coefs = coef_fxns(locs = sp_locs, coef = coef_comp)
    colnames(coefs) = c("alpha0","alpha1","b1","b2","b3","b4","b5")
    # Vector of the coefficient vectors beta(s1), beta(s2),..., beta(sn).
    B = coefs[,c("b1","b2","b3","b4","b5")]
    
    # Generating the true response data y
    for(i in 1:years){ # iterating over the ith obs at all locations
      scores = pc.scores_list[[i]] # the ith matrix of pc scores
      for(loc in 1:n_locs){ # iterating over spatial locations
        # One element of y at location loc and year i
        sim_nonfd[(loc-1)*years + i, "y"] = coefs[,c("alpha0")][loc] + coefs[loc,c("alpha1")]%*%sim_nonfd[(loc-1)*years + i,"z"] + t(scores[,loc])%*%B[loc,]
        sim_nonfd[(loc-1)*years + i, "CountyI"] = cnty_ids[loc]
        sim_nonfd[(loc-1)*years + i, c("long", "lat")] = sp_locs[loc,] # sp_locs[loc,] is a double
      }
    }
    
    # Calculating the noise for the response variable based on the signal to noise ratio
    if(is.null(sigma)){
      sigma = sqrt(var(sim_nonfd$y)/snr)
    }
    # y without noise
    y_noerr = sim_nonfd$y
    # Adding measurement error to y
    sim_nonfd$y = sim_nonfd$y + rnorm(n, sd = sigma)
    # # The tau-th quantile of the error distribution
    quant_eps = qnorm(tau, sd = sigma)
    # The true conditional quantile Q_Y(tau|X_i)
    true_cond_quant = y_noerr + quant_eps
    
    ############## Training and Testing for MSPE (K-fold CV) #################
    train_n.obs = ceiling((n - (n/n_folds))/n_locs) # number of training observations per county
    train.indx_mat = combn(1:years, train_n.obs)[,1:n_folds] 
    # Iterating over the number of folds
    folds = n_folds
    if (fld_lim < folds) folds = fld_lim
    for(m in 1:folds){
      cat("   fold", m, "\n")
      # The train and test indices for the mth fold
      train_indx = train.indx_mat[,m] 
      test_indx = c(1:years)[-train_indx]
      # A vector to store the testing ID's for this county 
      test_ids = c()
      # Iterating over the number of countys to create the training and testing sets
      for(cnty in unique(sim_nonfd$CountyI)){
        test_ids = append(test_ids, sim_nonfd[sim_nonfd$CountyI == cnty,][test_indx,]$Id)  # Add obs Id to tst_ids
      } # End of for loop for counties
  
      # The training and test data for nonfd and fundat using the Id's 
      test_nonfd = sim_nonfd[sim_nonfd$Id %in% test_ids,]
      test_fd = X_new[X_new$Id %in% test_ids,]
      train_nonfd = sim_nonfd[!sim_nonfd$Id %in% test_ids,]
      train_fd = X_new[!X_new$Id %in% test_ids,]
      smooth_trn_fd = fd_smooth_2steps_trn_tst(fd_dat_obj = train_fd, n_obs = nrow(train_fd)/N, n_days = N, sim = T, 
                                               reval = reval, llam.lims = log.lam)
      min_opt_lam = smooth_trn_fd$min_opt_lambda
      max_opt_lam = smooth_trn_fd$max_opt_lambda
      comb_train_fd_lst = smooth_trn_fd$comb_fd
      
      # Changing the column names of train_nonfd 
      colnames(train_nonfd) = colnames(sim_nonfd)
      
      # Using the optimal lambdas found in training set to smooth test set
      smooth_tst_fd = fd_smooth_2steps_trn_tst(test_fd, n_obs = nrow(test_fd)/N, min_opt_lamb = min_opt_lam, max_opt_lamb = max_opt_lam,
                                               n_days = N, sim = T, reval = reval, llam.lims = log.lam)
      comb_test_fd_lst = smooth_tst_fd$comb_fd
      colnames(test_nonfd) = colnames(sim_nonfd)
      
      # Splitting the true conditional quantiles into train and test
      train_cond_quant = true_cond_quant[!sim_nonfd$Id %in% test_ids]
      test_cond_quant = true_cond_quant[sim_nonfd$Id %in% test_ids]
      
      # The current index for the vectors/matrices
      run_indx = (n_folds*(iter-1)) + m
      ####### Fitting the svqfm model with the training data #########
      if("svqfm" %in% models){
        # Without penalty
        train_mod = train_fxn_qr(comb_train_fd_lst, train_nonfd, sim = T,
                                         nharm = nharm, spat_tri = spat_tri, d = d, cntr_fxns = T,
                                         true_quant = train_cond_quant, tau = tau, pnlty = pnlty, n_lam = n_lam, scl = scl,
                                         tmpmax_pars = tmpmax_pars, algo = qr_algo)
        svqfm_train_mse[run_indx] = train_mod$train_mse
        # Printing the train MSE
        cat("     The SVQFM train prediction error is", train_mod$train_mse, "\n")
        # Calculate the predicted quantiles using test data
        # With penalty 
        test_res = test_fxn_qr(comb_test_fd_lst, test_nonfd, train_mod$cnty_means, train_mod$fpca_obj,
                                      sim = T, num_harm = nharm, cntr_fxns = T, true_quant = test_cond_quant,
                                      spline_coefs = train_mod$spline_coefs, train_Q2 = train_mod$train_Q2,
                                      triang = spat_tri, tau = tau)
        # Returning the test MSE
        svqfm_mspe[run_indx] = test_res$test_MSE
        # Printing test MSE
        cat("     The SVQFM test prediction error is", svqfm_mspe[run_indx],"\n")
      }

      ############ Call the other function for FQR ##############
      if(comp_model){
        # Giving the training and testing data to fit the FQR
        comp_mod_res = model_comp_qr(train.lst_fd = comb_train_fd_lst, test.lst_fd = comb_test_fd_lst, nonfd.train = train_nonfd, 
                                     nonfd.test = test_nonfd, train_quant = train_cond_quant, test_quant = test_cond_quant, 
                                     num_harm = nharm, sim = T, tau_val = tau)
        if("fqr" %in% models){
          fqr_train_mse[run_indx] = comp_mod_res$train_mse
          fqr_mspe[run_indx] = comp_mod_res$test_mse
          cat("      The FQR train prediction error  is",comp_mod_res$train_mse, "\n")
          cat("      The FQR test prediction error  is",comp_mod_res$test_mse,"\n")
        }
      }
    } # end of for loop for folds
  } # End of for loop for repetitions
  
  # Storing the results for each of the models after all repetitions have been completed
  # List to store the results for all models
  model_res = list()
  # Spatially Varying Quantile Functional Model
  if("svqfm" %in% models){
    model_res[["svqfm"]] = list("train_mse" = svqfm_train_mse,"test_mse" = svqfm_mspe, "train_error_diff" = train_mod$error_diff, "test_error_diff" = test_res$error_diff,
                                "train_eta" = train_mod$train_etas)
  }
  # Functional Quantile Regression
  if("fqr" %in% models){
    model_res[["fqr"]] = list("train_mse" = fqr_train_mse,"test_mse" = fqr_mspe)
  }
  return(model_res)
}












