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
library(RColorBrewer)
library(ggplot2)
library(classInt)
### These packages are for handling the maps and location data ###
library(raster)
library(rnaturalearth)
library(sf)
library(rgdal)
library(rgeos)
library(tigris)
library(maps)

source_rmd = function(file, ...) {
  tmp_file = tempfile(fileext=".R")
  on.exit(unlink(tmp_file), add = TRUE)
  knitr::purl(file, output=tmp_file)
  source(file = tmp_file, ...)
}

# This function adds the county names to the functional data 
cnty_fxn = function(fdat, non_fd){
  # new column for county in functional data
  fdat[,"county"] = rep(NA, dim(fdat)[1])
  yrs = unique(fdat$Year)
  for(yr in yrs){
    reg_yr = non_fd[non_fd$Year==yr,]
    for(i in 1:length(reg_yr$CountyI)){
      id = reg_yr$CountyI[i]
      name = reg_yr$county[i]
      fdat$county[fdat$Year == yr & fdat$CountyI == id] = name 
    }
  }
  # After adding the county column, there are still rows with NA values
  # for county. This is because the number of counties in non_fd is not equivalent
  # to the number of counties in fun_dat even after removing the NA values.
  # Removing the NA values from fun_dat
  fdat_var = na.omit(fdat)
  # Changing the county names to lower case so that they match with non_fd.
  fdat_var$county = sapply(fdat_var$county, FUN = tolower)
  
  return(fdat_var)
}

# This function adds the county names to the functional data 
cnty_fxn_allfd = function(fdat, non_fd){
  # new column for county in functional data
  print(c("num obs in fdat Before removing temp NA", dim(fdat)[1]/365))
  fdat_nona_temp = na.omit(fdat)  # Remove all fdat with temp NA
  print(c("num obs in fdat After removing temp NA", dim(fdat_nona_temp)[1]/365))
  fdat_nona_temp[,"county"] = rep(NA, dim(fdat_nona_temp)[1])
  yrs = unique(fdat_nona_temp$Year)
  # Adding the county name for each county ID
  for(cid in unique(non_fd$CountyI)){
    name = non_fd[non_fd$CountyI == cid,"county"][1]
    yrs = non_fd[non_fd$CountyI == cid, "Year"]
    fdat_nona_temp[fdat_nona_temp$CountyI == cid & fdat_nona_temp$Year %in% yrs,"county"] = name
  }
  # Changing the county names to lower case so that they match with non_fd.
  fdat_with_cnty = fdat_nona_temp[!is.na(fdat_nona_temp$county),]
  fdat_without_cnty = fdat_nona_temp[is.na(fdat_nona_temp$county),]
  # Stack rows such that all fundat that has county (and are in regdat) will be first, followed by fundat w/o county
  fdat_final <- rbind(fdat_with_cnty, fdat_without_cnty)

  return(fdat_final)
}

center_fxn = function(st, st_bnd, state_centroid){
  # st: states
  # st_bnd: the boundaries of the states
  # state_centroid: the center of the boundary
  
  # list of the centroids and counties for each state
  center_list = list()
  st_cntys_list = list()
  for(i in 1:length(st)){
    st_poly = counties(st[i], cb = TRUE) # the multipolygon for ith state
    st_cntys_list[[i]] = c(tolower(st_poly$NAME)) # the counties corresponding to the centroids
    # Find the centroids
    centroids = st_centroid(st_geometry(st_poly))
    # Creating a matrix for the counties centroids
    num_cnt = length(st_poly$NAME)
    centr_mat = matrix(0, nrow = num_cnt, ncol = 2)
    for(j in 1:num_cnt){
      centr_mat[j,] = centroids[[j]]
    }
    # Centering the data based on the centroid for the given boundary
    center_list[[i]] = centr_mat - matrix(state_centroid, nrow = dim(centr_mat)[1], ncol = 2, byrow = T) # centering the coordinates
  }
  return(list(center_list, st_cntys_list))
}

fd_smooth = function(fd_dat_obj, n_obs, n_knots = 25, n_days = 365,
                     reval = c(0,1), llam.lims = c(-15, 5), lam.len_ = 50, 
                     dfscale_ = 1){
  min_temp = matrix(fd_dat_obj$TMIN, nrow = n_obs, ncol = n_days, byrow = TRUE)
  max_temp = matrix(fd_dat_obj$TMAX, nrow = n_obs, ncol = n_days, byrow = TRUE)
  # matrix of time points for all observations
  time_pts = matrix(seq(0, 1, len = n_days), nr = n_obs, nc = n_days, byrow = T) 
  # mintrn.fd and maxtrn.fd contain the list of smoothed data, and other parameters for X1 and X2
  fd1 = dat2fd2(argvals = t(time_pts), y = t(min_temp), rangeval = reval, bs = NULL,
                loglam.lims = llam.lims, lam.len = lam.len_, dfscale = dfscale_, nk = n_knots)
  fd2 = dat2fd2(argvals = t(time_pts), y = t(max_temp), rangeval = reval, bs = NULL,
                loglam.lims = llam.lims, lam.len = lam.len_, dfscale = dfscale_, nk = n_knots)
  # The smoothed functional data for X1 and X2
  fd_X1 = fd1$fd # for the min temp
  fd_X2 = fd2$fd # for the max temp
  # The combined X1 and X2 data into a multivariate data object
  fd_comb = combine.fd(fd_X1, fd_X2) 
  return(fd_comb)
}

fd_smooth_2steps_trn_tst = function(fd_dat_obj, n_obs, min_opt_lamb = NULL,
                                    max_opt_lamb = NULL, n_knots = 50, n_days = 365,
                                    reval = c(1,365), llam.lims = c(-15.5,-8), lam.len_ = 50,
                                    dfscale_ = 1, sim = F){
  # The unique counties
  uniq_cntId = unique(fd_dat_obj$CountyI)
  # First take ALL samples from ALL locations and ALL years and find the best lambda for smoothing
  print_exld = F
  if (is.null(min_opt_lamb)) {
    print_exld = T
    # If this is not a simulation, the data should be organized into a matrix.
    time_pts = matrix(seq(reval[1],reval[2], len = n_days), nr = n_obs, nc = n_days, byrow = T)

    if(!sim){
      min_temp = matrix(fd_dat_obj$TMIN, nrow = n_obs, ncol = n_days, byrow = TRUE)
      max_temp = matrix(fd_dat_obj$TMAX, nrow = n_obs, ncol = n_days, byrow = TRUE)
    } else{
      # If it's a simulation, the variables are called x1 and x2 for min and max temperature
      min_temp = matrix(fd_dat_obj$x1, nrow = n_obs, ncol = n_days, byrow = TRUE)
      max_temp = matrix(fd_dat_obj$x2, nrow = n_obs, ncol = n_days, byrow = TRUE)
    }
    # mintrn.fd and maxtrn.fd contain the list of smoothed data, and other parameters for X1 and X2
    min_opt_lamb = dat2fd2_opt_lambda(argvals = t(time_pts), y = t(min_temp), rangeval = reval, bs = NULL,
                                      loglam.lims = llam.lims, lam.len = lam.len_, dfscale = dfscale_, nk = n_knots)
    max_opt_lamb = dat2fd2_opt_lambda(argvals = t(time_pts), y = t(max_temp), rangeval = reval, bs = NULL,
                                      loglam.lims = llam.lims, lam.len = lam.len_, dfscale = dfscale_, nk = n_knots)
  }
  # Now that we have best lambda for min and for max we use them to smooth sample on each location
  # and putting the smoothed data by county into a list so we can find mean function per location.
 
  # Smoothing the temp data for each county
  fds_comb = list()
  for (i in 1:length(uniq_cntId)){
    cntyId = uniq_cntId[i]
    num_obs = nrow(fd_dat_obj[fd_dat_obj$CountyI == cntyId,])/n_days
    time_pts_cnty = matrix(seq(reval[1],reval[2], len = n_days), nr = num_obs, nc = n_days, byrow = T)
    # Adding in an if condition to account for whether or not this is a simulation
    if(!sim){
      # Find the indices all fd_dat_obj that are in CountyI but has county NA
      countyI_with_cnty_na_inds = is.na(fd_dat_obj[fd_dat_obj$CountyI == cntyId,]$county)
      # Create a year indices for which county is NA, these are years for the county
      # for which there is no nonfd data
      countI_yr_with_cnty_na_inds = which(countyI_with_cnty_na_inds[seq(1, length(countyI_with_cnty_na_inds), 365)])
      fd_min_cnty = fd_dat_obj[fd_dat_obj$CountyI == cntyId,]$TMIN
      fd_max_cnty = fd_dat_obj[fd_dat_obj$CountyI == cntyId,]$TMAX
      min_temp_cnty = matrix(fd_min_cnty, nrow = num_obs, ncol = n_days, byrow = TRUE)
      max_temp_cnty = matrix(fd_max_cnty, nrow = num_obs, ncol = n_days, byrow = TRUE)
    } else{
      fd_min_cnty = fd_dat_obj[fd_dat_obj$CountyI == cntyId,]$x1
      fd_max_cnty = fd_dat_obj[fd_dat_obj$CountyI == cntyId,]$x2
      min_temp_cnty = matrix(fd_min_cnty, nrow = num_obs, ncol = n_days, byrow = TRUE)
      max_temp_cnty = matrix(fd_max_cnty, nrow = num_obs, ncol = n_days, byrow = TRUE)
    }
    
    fd1 = dat2fd2_withLambda(argvals = t(time_pts_cnty), y = t(min_temp_cnty), rangeval = reval, bs = NULL,
                             lambda = min_opt_lamb, dfscale = dfscale_, nk = n_knots)
    fd2 = dat2fd2_withLambda(argvals = t(time_pts_cnty), y = t(max_temp_cnty), rangeval = reval, bs = NULL,
                             lambda = max_opt_lamb, dfscale = dfscale_, nk = n_knots)
    # The smoothed functional data for X1 and X2
    fd_X1 = fd1 # for the min temp
    fd_X2 = fd2 # for the max temp
    # The combined X1 and X2 data into a multivariate data object
    fd_comb = combine.fd(fd_X1, fd_X2)
    fd_comb$CountyI = cntyId   # Added county ID
    if(!sim) fd_comb$excld_yrs_frm_trn = countI_yr_with_cnty_na_inds
    fds_comb[[i]] = fd_comb
  }
  return(list("comb_fd" = fds_comb, "min_opt_lambda" = min_opt_lamb, "max_opt_lambda" = max_opt_lamb))
}



# Function that performs FPCA
fpca = function(fun_dat_obj, trn, num_harm){
  # The training data to perform FPCA
  fd_train = fun_dat_obj[trn]
  # Joint FPCA of the smoothed functional data X1 and X2
  fpca_train = joint.fpca(fd_train, nharm = num_harm)
  # The joint functional principal component scores
  fpc.scores_train = fpca_train$pc_scores
  # I need to return the fpca object so that I can calculate the scores for 
  # the training data
  return(list("scores" = fpc.scores.train, "fpca_obj" = fpca_train$fpca))
}

tri_fxn = function(sp_bnd, n_tri){
  # sp_bnd: A matrix containing the coordinates for the spatial boundary.
  # n_tri: An integer that controls the fineness of the triangulation.
  
  # Creates a triangulation of the entire spatial domain
  return(tri_mat = TriMesh(Pt = sp_bnd, n = n_tri))
}

dat2fd2_opt_lambda <- function(argvals, y, rangeval, bs=NULL, loglam.lims=c(-15,5), lam.len=50,
                               dfscale=1, traceit=F, norder=4, nk=50){
  #############################################################
  # only return the best lambda found by computing gcvs for all
  # data points over all location all years, among a sequence of lambda.
  #############################################################
  
  # argvals: the time points for all functional observations
  # y: the simulated functional data in the form of a matrix
  # rangeeval : lower and upper end of time interval
  # bs: this is just a placeholder variable for an if statement
  # loglam.lims: this the interval for the smoothing parameter
  # lam.len: the number of smoothing parameters to test.
  # dfscale: not sure what this is, doesn't show up in the function
  # traceit: also some placeholder variable so an if statement can work
  # norder: the order of the B-spline basis
  # nk: the number of knots for the b-spline basis.

  if (is.null(bs)){
    norder <- norder
    nk <- nk
    #x <- sort(argvals); n <- length(x); nk <- 20
    #times <- quantile(x[2:(n-1)], seq(0, 1, length=nk)) # knots
    # this is the knots for the b-spline
    times <- seq(rangeval[1], rangeval[2], length=nk)
    # This is the way they calculate the num of basis functions if only the number
    # of breaks is given.
    nbasis <- length(times) + norder - 2
    bs <- create.bspline.basis(rangeval, nbasis, norder, breaks=times)
  }

  ###
  if (nrow(unique(t(argvals)))>1){
    stop("the smooth.basis from fda package does not perform correctly if time points are different")
  } else {
    argvals <- argvals[,1]
  }

  loglam <- seq(loglam.lims[1], loglam.lims[2], len=lam.len)
  gcvs <- array(dim=lam.len)
  mmse <- array(dim=lam.len)
  dfs <- array(dim=lam.len)
  sobjs <- list()
  for (i in (1:lam.len)){
    if (traceit) cat(i, "\n")
    lambdai <- exp(loglam[i])
    fP <- fdPar(bs, 2, lambdai)
    sobjs[[i]] <- smooth.basis(argvals=argvals, y=y, fdParobj=fP)
    gcvs[i] <- mean(sobjs[[i]]$gcv)
    mmse[i] <- sobjs[[i]]$SSE/prod(dim(y))
    dfs[i] <- mean(sobjs[[i]]$df)
  }
  # I think here he is finding which smoothing parameter gives the smallest gcv 
  # and then the corresponding smooth.basis object is sobj.
  gcvId <- which.min(gcvs)
  
  return(lambda=exp(loglam[gcvId]))
}

dat2fd2_withLambda <- function(argvals, y, rangeval, bs=NULL, lambda,
                               dfscale=1, traceit=F, norder=4, nk=50){
  #############################################################
  # Compute smoothed data given the optimal lambda found earlier 
  # using the dat2fd2_opt_lambda function.
  #############################################################
  
  # argvals: the time points for all functional observations
  # y: the simulated functional data in the form of a matrix
  # rangeeval : lower and upper end of time interval
  # bs: this is just a placeholder variable for an if statement
  # lambda: the lambda with which to perform smoothing
  # dfscale: not sure what this is, doesn't show up in the function
  # traceit: also some placeholder variable so an if statement can work
  # norder: the order of the B-spline basis
  # nk: the number of knots for the b-spline basis.
  if (is.null(bs)){
    norder <- norder
    nk <- nk
    # this is the knots for the b-spline
    times <- seq(rangeval[1], rangeval[2], length=nk)
    # This is the way they calculate the num of basis functions if only the number
    # of breaks is given.
    nbasis <- length(times) + norder - 2
    bs <- create.bspline.basis(rangeval, nbasis, norder, breaks=times)
  }
  ###
  if (nrow(unique(t(argvals)))>1){
    stop("the smooth.basis from fda package does not perform correctly if time points are different")
  } else {
    argvals <- argvals[,1]
  }

  fP <- fdPar(bs, 2, lambda)
  sobj <- smooth.basis(argvals=argvals, y=y, fdParobj=fP)
  
  return(fd=sobj$fd)
}

# I copied this function from the file "plfam.R". It performs 
# This function combines the two functional data objects into one 
# multivariate data object using an array.
combine.fd <- function(fdobj1, fdobj2){
  bscoef <- array(dim=c(dim(fdobj1$coefs),2))
  bscoef[,,1] <- fdobj1$coefs
  bscoef[,,2] <- fdobj2$coefs
  fdobj <- fd(bscoef, fdobj1$basis)
  return(fdobj)
}

# This function performs multivariate PCA.
joint.fpca <- function(fdobj, nharm = NULL, thresh = 0.80, centerfns = T, simulation=F, scl = T){
  # Setting the number of harmonics to the number of harmonics for the smoothing splines
  # of the fd data. I
  fpcaobj <- pca.fd(fdobj, nharm = fdobj$basis$nbasis, centerfns=centerfns)
  # Picking the number of harmonics based on the percentage of variance that they
  # explain.
  if (is.null(nharm)) nharm <- sum(cumsum(fpcaobj$varprop) <= thresh)
  if (length(dim(fpcaobj$scores)) <= 2){
    score <- fpcaobj$scores
  } else {
    # When you do joint FPCA, you sum the PC scores for the the two 
    # variables. Scores is an array so we are summing with respect to 
    # the row and column (this is c(1,2)). So we end up with only one 
    # matrix of PC scores as opposed to 2 (or however many variables there are).
    # Refer to equation (8.21) in Ramsey FDA book to understand what they're doing here.
    score <- apply(fpcaobj$scores, c(1,2), sum)
  }
  if (!simulation) cat("The number of harmonics used for", thresh, " variations is", nharm, "\n")
  return(list(fpca = fpcaobj, pc_scores = score[, 1:nharm]))
}

# For a new obs, we project the data onto the eigenfunctions calculated
# using the training data.
proj.mfpca_nv <- function(pcafd, fdobj, cntr = F){
  # pcafd - the fpc obj from training
  # fdobj - the new functional data object in Bspline basis form and after X_1 and X_2 combined
  if (!(pcafd$harmonics$basis == fdobj$basis)){
    stop("basis mismatched.")
  }
  
  coef <- fdobj$coefs  # spline coefs of new data
  coefd <- dim(coef)
  nvar <- coefd[3]
  nrep <- coefd[2]    # I think this is the number of functional vars
  nharm <- dim(pcafd$harmonics$coefs)[2]  
  
  # Checking if the functional data should be centered or not
  if(cntr){
    # de-mean by the pcafd$meanfd
    for (i in (1:nvar)){    # Removing the training data mean coeffs from each new obs coeffs
      fdobj$coefs[,,i] <- fdobj$coefs[,,i] - pcafd$meanfd$coefs[,,i]
    }
  }
  
  harmscr <- array(0, c(nrep, nharm, nvar))
  coefarray <- fdobj$coefs  # spline coefs of new data
  harmcoefarray <- pcafd$harmonics$coefs   # spline coeffs of eigen function from training data
  basisobj <- fdobj$basis   # Spline basis common to training and new data
  for (j in 1:nvar){
    fdobjj <- fd(as.matrix(coefarray[, , j]), basisobj)  # fd obj of newdata
    harmfdj <- fd(as.matrix(harmcoefarray[, , j]), basisobj)  # fd obj of eigen functions
    harmscr[, , j] <- inprod(fdobjj, harmfdj)   # inner product which is projection of obs on each eigen
  }
  score <- apply(harmscr, c(1,2), sum)

  return (list(harmscr = harmscr, pc_scores = score))
}

preprocessing_nosmooth = function(fun_dat, non_fd, states, vars, st_bnd,
                         n_knts = 25, n_days = 365, reval = c(0,1), llam.lims = c(-15, 5), 
                         lam.len_ = 50, dfscale_ = 1){
  # The arguments for this function are
  # fun_dat: the functional data object
  # non_fd: the scalar predictors (including response)
  # states: list of states
  # vars: the variables of interest 
  # n_tri: controls the fineness of the triangulation
  # st_bnd: a matrix where each row corresponds to each states boundaries.
  # trn_split: the split for the training and testing
  # There are many arguments which are simply smoothing parameters for the functional
  # data. 
  
  # Changing the column name for county in the non_fd data
  # need to include an if statement here to change the column name "County" to county
  for(i in 1:length(colnames(non_fd))){
    if(colnames(non_fd)[i] == "County") colnames(non_fd)[i] = "county"
  }
  # Changing lower case for non_fd county
  non_fd$county = tolower(non_fd$county)
  # Changing the county in non_fd to lower case for comparison later on
  non_fd$State = tolower(non_fd$State)
  
  ############# Adding the county names to fun_dat ############# 
  fd_dat = cnty_fxn_allfd(fun_dat, non_fd)
  # fd_dat contain data for county/year that are in regdat as well as not in regdat
  # The reason to include those that are not is to use more data for FPCA
  
  ########## Finding the centroids of all states counties #####################
  # Before calling the center_fxn, calculate the center of the given boundary, st_bnd
  # and then give it as an argument to center_fxn
  st_mat = st_bnd[[1]][,c("long","lat")] # initializing the matrix
  # First need to bind the matrices in the list st_bnd
  if(length(st_bnd) > 1){
    for(s in 2:length(st_bnd)){
      st_mat = rbind(st_mat, st_bnd[[s]][,c("long","lat")])
    }
  }
  # Calculating the center of the boundary (the entire Midwest boundary)
  st_lon_cent = (min(st_mat[,"long"]) + max(st_mat[,"long"]))/2
  st_lat_cent = (min(st_mat[,"lat"]) + max(st_mat[,"lat"]))/2
  st_cent = c(st_lon_cent, st_lat_cent)
  # vector of the counties and centroids
  cntr_fxn_otpt = center_fxn(states, st_bnd, state_centroid = st_cent)
  cntrs = cntr_fxn_otpt[[1]] # The list of centroid matrices for each state
  cntys = cntr_fxn_otpt[[2]] # The list of the county names
  # Removing spaces in county name
  non_fd$county = gsub(" ", "", non_fd$county)
  fd_dat$county = gsub(" ", "", fd_dat$county)
  # For each state, I am adding the corresponding county centers
  for(l in 1:length(states)){
    # For each states county names (cntys[[l]]), remove all punctuation and white space from the county name so that I can 
    # compare it to the county names in fd and nonfd. 
    cnty_vec = gsub(" |'|\\.", "", cntys[[l]])
    for(k in 1:nrow(cntrs[[l]])){
      # For each county, find the corresponding county Id using the state and county name from cntys[[l]][k]
      # Adding the centroids to functional data and the non_fd data
      # Changed states[l] to lower case so that it is comparable with non_fd$State
      cnty_id = non_fd[non_fd$State == tolower(states[l]) & non_fd$county == cnty_vec[k],"CountyI"][1]
      if (is.na(cnty_id)) next   # No need to update centroid since we don't have this county in this state
      # The center for the kth county in the lth state
      cntr = data.frame(matrix(cntrs[[l]][k,], nrow = 1, ncol = 2))
      # Include a condition to check the county Id since some states have the same county names
      fd_dat[!is.na(fd_dat$county) & fd_dat$county == cnty_vec[k] & fd_dat$CountyI == cnty_id, c("long","lat")] = cntr
      non_fd[non_fd$county == cnty_vec[k] & non_fd$CountyI == cnty_id, c("long","lat")] = cntr
    }
  }
  
  ################## Organizing the scalar predictors and functional data ###############
  mod_vars = c(vars, "long","lat", "county")
  pred_vars = data.frame(matrix(0, nrow = dim(non_fd)[1], ncol = length(mod_vars)))
  i = 1
  j = 1
  # The while loop iterates over the total number of replicates in the non_fd dataset.
  # For the jth replicate, we want to compare the year and county in non_fd to the 
  # ith replicate's year and county in fun_dat_var. For each replicate in fun_dat_var, the 
  # year and county is the same so we only need to iterate every 365, the number 
  # of time points for each replicate. For each comparison, we add the the corresponding
  # predictor variables from non_fd. 
  while(j < nrow(pred_vars) + 1){
    # need to add a county here for the county Id since there are some states 
    # that have the same county names
    row_ind = which(non_fd$Year == fd_dat[i,]$Year & !is.na(fd_dat[i,]$county)
                    & non_fd$county == fd_dat[i,]$county & non_fd$CountyI == fd_dat[i,]$CountyI)
    if(is.na(fd_dat[i,]$county)){
      print(row_ind)
      print(non_fd[row_ind, mod_vars])
    }
    pred_vars[j,] = non_fd[row_ind, mod_vars]
    i = i + 365
    j = j + 1
  }
  
  # Organizing the variables
  if("Irrigated.eprop" %in% mod_vars & length(states) == 1) colnames(pred_vars) = c(mod_vars[1:2], "irrig_prop", mod_vars[4:length(mod_vars)])
  # Adding a condition here for the Midwest since the variable name for irrigation proportion is different
  else if ("Irrigated.actualprop" %in% mod_vars) colnames(pred_vars) = c(mod_vars[1:2], "irrig_prop", mod_vars[4:length(mod_vars)])
  else colnames(pred_vars) = mod_vars

  return(list("non_fd" = pred_vars, "fd_dat" = fd_dat))
}

