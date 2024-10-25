library(MASS)
library(dplyr)
library(tidyr)
source("GSVFM/util_funcs/plfam_yehua.R")
source("GSVFM/PLFAM_FLM_ExtractID_Midwest.R")

################################################
#### additional functions for data analysis ####
################################################
getId <- function(id1, id2){
  # return the order of id2 to match id1
  out <- array(dim=length(id1))
  for (i in (1:length(id1))){
    temp <- which(id2==id1[i])
    if (length(temp)>1){
      stop("more than one id2 element correspond to id[i]")
    } else {
      out[i] <- which(id2==id1[i])
    }
  }
  return(out)
}

get.sobj <- function(fundat, ii, loglam.lims= c(3,5), lam.len=50, nk=50, reval = c(0,1)){ # ii: 1=TMAX 2=TMIN
  if (ii==1){
    temp <-spread(fundat[, c("Year", "DateI", "Id", "TMAX")], DateI, TMAX)
  } else if (ii==2){
    temp <-spread(fundat[, c("Year", "DateI", "Id", "TMIN")], DateI, TMIN)
  }
  allX <- t(temp)[-c(1,2),] # contain NA
  X <- allX[,!is.na(allX[1,])] # remove NA
  id <- temp$Id[!is.na(allX[1,])]
  year <- temp$Year[!is.na(allX[1,])]
  argvals <- matrix(seq(reval[1],reval[2],len = 365), nr=365, nc=ncol(X))
  sobj <- dat2fd2(argvals=argvals, y=X, rangeval=reval,
                  loglam.lims=loglam.lims, lam.len=lam.len, nk=nk)
  return(list(id=id, year=year, sobj=sobj))
}

demean <- function(fd, year){
  y <- unique(year)
  mfd <- list()
  for (i in (1:length(y))){
    ind <- year==y[i]
    mfd[[i]] <- mean(fd[ind])
    fd$coefs[,ind] <- fd$coefs[,ind] - as.vector(mfd[[i]]$coefs)
  }
  return(list(mfd=mfd, fd=fd))
}

get.sobj2 <- function(fundat, ii, lam, basis, reval = c(0,1)){ # ii: 1=TMAX 2=TMIN
  if (ii==1){
    temp <-spread(fundat[, c("Year", "DateI", "Id", "TMAX")], DateI, TMAX)
  } else if (ii==2){
    temp <-spread(fundat[, c("Year", "DateI", "Id", "TMIN")], DateI, TMIN)
  }
  allX <- t(temp)[-c(1,2),] # contain NA
  X <- allX[,!is.na(allX[1,])] # remove NA
  id <- temp$Id[!is.na(allX[1,])]
  year <- temp$Year[!is.na(allX[1,])]
  argvals <- matrix(seq(reval[1],reval[2], len = 365), nr=365, nc=ncol(X))
  fP <- fdPar(basis, 2, lam)
  sobj <- smooth.basis(argvals=argvals, y=X, fdParobj=fP)
  return(list(id=id, year=year, sobj=sobj))
}

####################
#### containers ####
####################
##### Arrays to store the PLFAM prediction errors #######
pe2all <- array(dim=c(n_iters, n_folds)) # PE of PLFAM(joint)
wpe2all <- array(dim=c(n_iters, n_folds)) # WPE of PLFAM(joint)
pe2all_mse_perc <- array(dim=c(n_iters, n_folds))
wpe2all_mse_perc  <- array(dim=c(n_iters, n_folds))
pe2all_mse_var <- array(dim=c(n_iters, n_folds))
wpe2all_mse_var <- array(dim=c(n_iters, n_folds))
pe2all_mse_mean_abs <- array(dim=c(n_iters, n_folds))
wpe2all_mse_mean_abs <- array(dim=c(n_iters, n_folds))
pe2all_mse_mean_sqr <- array(dim=c(n_iters, n_folds))
wpe2all_mse_mean_sqr <- array(dim=c(n_iters, n_folds))
plfam_train_time <- array(dim=c(n_iters, n_folds))

##### Arrays to store the FLM prediction errors #######
pe2linearall <- array(dim=c(n_iters, n_folds)) # PE of FLM-Cov(joint)
wpe2linearall <- array(dim=c(n_iters, n_folds)) # WPE of FLM-Cov(joint)
pe2linearall_mse_perc <- array(dim=c(n_iters, n_folds))
wpe2linearall_mse_perc <- array(dim=c(n_iters, n_folds))
pe2linearall_mse_var <- array(dim=c(n_iters, n_folds))
wpe2linearall_mse_var <- array(dim=c(n_iters, n_folds))
pe2linearall_mse_mean_abs <- array(dim=c(n_iters, n_folds))
wpe2linearall_mse_mean_abs <- array(dim=c(n_iters, n_folds))
pe2linearall_mse_mean_sqr <- array(dim=c(n_iters, n_folds))
wpe2linearall_mse_mean_sqr <- array(dim=c(n_iters, n_folds))
totalweight2 <- array(dim=c(n_iters, n_folds))
flm_train_time <- array(dim=c(n_iters, n_folds))

# Data using pre-processing function
load("Data_Extraction_cleaning/RdataFiles/RawMidwestFundatPreProc_incl_irrig.RData")
load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcRegdat.RData")
# The midwest data 
fundat0 = midw_fd_new
regdat0 = midw_regdat_new

n_iters = 9 # number of repetitions
n_folds <- 5 # number of folds in validations
DE_MEAN_RESP = T
llam.lims = c(-15.5,-8)

#############################
#### code for one period ####
#############################

# define expression for tryCatch
oneiteration <- expression({
  set.seed(12876+i*123)
  for(vv in 1:n_folds){
    cat("\n","Fold", vv, "\n")
    cat('------------', "\n")
    # There are more rows in test_inds than there are test observations
    # just because I generated the matrix with 3000 rows.
    Id.valid <- na.omit(test_inds[,(i-1)*n_folds + vv])
    fundat1 <- filter(fundat0, !(Id %in% Id.valid)) # note it may contain more counties than trainID
    # Including the code with excl_cnty_ids since when I choose the train indices, I take everything that's not in the test 
    # indices which can include countys that have less than 5 years of data. 
    fundat1 <- filter(fundat1, !fundat1$CountyI %in% excl_cnty_ids) # train fd 
    fundat2 <- filter(fundat0, (Id %in% Id.valid))                  # test fd
    regdat1 <- filter(regdat0, !(Id %in% Id.valid))                 # train nonfd
    regdat1 <- filter(regdat1, !regdat1$CountyI %in% excl_cnty_ids) 
    regdat2 <- filter(regdat0, (Id %in% Id.valid))                  # test nonfd
    ### Centering the yield data by the mean per year ###
    # Add column yld_mean_per_yr to regdat1
    regdat1$yld_mean_per_yr = ave(regdat1$Yield, regdat1$Year)
    # add De-Meaned column to regdat1
    regdat1$de_meaned_Yield = regdat1$Yield - regdat1$yld_mean_per_yr
    # Add yld_mean_per_yr to regdat2
    regdat2$yld_mean_per_yr = NA
    for (yr in unique(regdat2$Year)){
      regdat2[regdat2$Year == yr,]$yld_mean_per_yr = unique(regdat1[regdat1$Year == yr,]$yld_mean_per_yr)
    }
    # Add column for De-meaned Yield in regdat2
    regdat2$de_meaned_Yield = regdat2$Yield - regdat2$yld_mean_per_yr
    
    #######################
    #### pre-smoothing ####
    #######################
    # functional predictors: TMAX, TMIN (separately)
    # using all years
    ## smooth using the same smoothing parameters
    #### TMAX
    temp <- get.sobj(fundat1, ii=1, reval = c(0,1), loglam.lims = llam.lims)
    id.TMAX <- temp$id
    year.TMAX <- temp$year
    sobj.TMAX <- temp$sobj
    sobj.TMAX$lambda <- temp$sobj$lambda
    #### TMIN
    temp <- get.sobj(fundat1, ii=2, reval = c(0,1), loglam.lims = llam.lims)
    id.TMIN <- temp$id
    year.TMIN <- temp$year
    sobj.TMIN <- temp$sobj
    sobj.TMIN$lambda <- temp$sobj$lambda
    
    ########################
    #### FPCA and PLFAM ####
    ########################
    if(!DE_MEAN_RESP) y1 <- regdat1$Yield
    else y1 = regdat1$de_meaned_Yield
    #### Including interaction term between irrigation and prcp ####
    U1 <- cbind(regdat1$irrig_prop, regdat1$avgPRCP,
                regdat1$irrig_prop*regdat1$avgPRCP)
    weight1 <- regdat1$Area
    
    #### joint #####
    # check id.TMIN and id.TMAX
    if (!any(id.TMIN!=id.TMAX)){
      # both sobj.TMAX and sobj.TMIN are constructed on the same basis
      fd.joint <- combine.fd(sobj.TMAX$fd, sobj.TMIN$fd) # because id.TMIN==id.TMAX
      t_fpca_start = Sys.time()
      zobj.joint <- get.zeta(fd.joint)
      t_fpca = Sys.time()-t_fpca_start 
      xi.joint1 <- zobj.joint$zeta[getId(regdat1$Id, id.TMAX),] # because id.TMIN==id.TMAX
      t_plfam_start = Sys.time()
      ##### PLFAM Fitting ########
      fit2 <- plcosso(y1, U1, xi.joint1, weight=weight1/sum(weight1), traceit=T)
      t_plfam_fit = Sys.time() - t_plfam_start
      #### FLM Fitting #######
      t_flm_start = Sys.time()
      fit2.linear <- myridge(y1, U1, xi.joint1, weight=weight1/sum(weight1))
      t_flm_fit = Sys.time() - t_flm_start
    }
    # Total training run time for each fold for plfam and flm
    t_plfam <- t_plfam_fit + t_fpca
    t_flm <- t_flm_fit + t_fpca
    
    ##########################
    #### prediction error ####
    ##########################
    if(!DE_MEAN_RESP) y2 <- regdat2$Yield
    else y2 = regdat2$de_meaned_Yield
    #### Including interaction term between irrigation and prcp ####
    U2 <- cbind(regdat2$irrig_prop, regdat2$avgPRCP,
                regdat2$irrig_prop*regdat2$avgPRCP)
    weight2 <- regdat2$Area
    
    temp <- get.sobj2(fundat2, ii=1, lam=sobj.TMAX$lambda, basis=sobj.TMAX$fd$basis,reval = c(0,1))
    id.TMAX2 <- temp$id
    sobj.TMAX2 <- temp$sobj
    
    temp <- get.sobj2(fundat2, ii=2, lam=sobj.TMIN$lambda, basis=sobj.TMIN$fd$basis,reval = c(0,1))
    id.TMIN2 <- temp$id
    sobj.TMIN2 <- temp$sobj
    
    #### joint #####
    # both sobj.TMAX and sobj.TMIN are constructed on the same basis
    fd.joint2 <- combine.fd(sobj.TMAX2$fd, sobj.TMIN2$fd) # because id.TMIN==id.TMAX
    zobj.joint2 <- get.zeta2(fd.joint2, fd.joint)
    xi.joint2 <- zobj.joint2$zeta[getId(regdat2$Id, id.TMAX2),] # because id.TMIN==id.TMAX
    
    ##### PLFAM Prediction Error ####
    pred2 <- predict.plcosso(fit2, U2, xi.joint2)
    pe2 <- sum((y2-pred2)^2)
    wpe2 <- sum(weight2*(y2-pred2)^2)
    pe2_mse_perc = sum(((y2-pred2)^2)/(y2^2))
    wpe2_mse_perc = sum(weight2*((y2-pred2)^2)/(y2^2))
    pe2_mse_var = sum(((y2-pred2)^2)/var(y2))
    wpe2_mse_var = sum(weight2*((y2-pred2)^2)/var(y2))
    pe2_mse_mean_abs = sum(((y2-pred2)^2)/mean(abs(y2)))
    wpe2_mse_mean_abs = sum(weight2*((y2-pred2)^2)/mean(abs(y2)))
    pe2_mse_mean_sqr = sum(((y2-pred2)^2)/mean(y2^2))
    wpe2_mse_mean_sqr = sum(weight2*((y2-pred2)^2)/mean(y2^2))

    ######## FLM prediction error #######
    pred2 <- as.vector(cbind(1, U2, xi.joint2) %*% fit2.linear$beta)
    pe2.linear <- sum((y2-pred2)^2)
    wpe2.linear <- sum(weight2*(y2-pred2)^2)
    pe2.linear_mse_perc = sum(((y2-pred2)^2)/(y2^2))
    wpe2.linear_mse_perc = sum(weight2*((y2-pred2)^2)/(y2^2))
    pe2.linear_mse_var = sum(((y2-pred2)^2)/var(y2))
    wpe2.linear_mse_var = sum(weight2*((y2-pred2)^2)/var(y2))
    pe2.linear_mse_mean_abs = sum(((y2-pred2)^2)/mean(abs(y2)))
    wpe2.linear_mse_mean_abs = sum(weight2*((y2-pred2)^2)/mean(abs(y2)))
    pe2.linear_mse_mean_sqr = sum(((y2-pred2)^2)/mean(y2^2))
    wpe2.linear_mse_mean_sqr = sum(weight2*((y2-pred2)^2)/mean(y2^2))

    ## PLFAM ##
    pe2all[i, vv] <- pe2
    wpe2all[i, vv] <- wpe2
    pe2all_mse_perc[i, vv] <- pe2_mse_perc
    wpe2all_mse_perc[i, vv] <- wpe2_mse_perc
    pe2all_mse_var[i, vv] <- pe2_mse_var
    wpe2all_mse_var[i, vv] <- wpe2_mse_var
    pe2all_mse_mean_abs[i, vv] <- pe2_mse_mean_abs
    wpe2all_mse_mean_abs[i, vv] <- wpe2_mse_mean_abs
    pe2all_mse_mean_sqr[i, vv] <- pe2_mse_mean_sqr
    wpe2all_mse_mean_sqr[i, vv] <- wpe2_mse_mean_sqr
    # PLFAM train time
    plfam_train_time[i,vv] <- t_plfam
    
    ## FLM ##
    pe2linearall[i, vv] <- pe2.linear
    wpe2linearall[i, vv] <- wpe2.linear
    pe2linearall_mse_perc[i, vv] <- pe2.linear_mse_perc
    wpe2linearall_mse_perc[i, vv] <- wpe2.linear_mse_perc
    pe2linearall_mse_var[i, vv] <- pe2.linear_mse_var
    wpe2linearall_mse_var[i, vv] <- wpe2.linear_mse_var
    pe2linearall_mse_mean_abs[i, vv] <- pe2.linear_mse_mean_abs
    wpe2linearall_mse_mean_abs[i, vv] <- wpe2.linear_mse_mean_abs
    pe2linearall_mse_mean_sqr[i, vv] <- pe2.linear_mse_mean_sqr
    wpe2linearall_mse_mean_sqr[i, vv] <- wpe2.linear_mse_mean_sqr
    # FLM train time
    flm_train_time[i, vv] <- t_flm
      
    totalweight2[i, vv] <- sum(weight2)
}})

############################
#### run the experiment ####
############################

##########################################################################################
# The model fitting above was originally split into two files where one ran iterations
# 1 to 4 and the other iterations 5 to 9. The code is the same but the the for loop over 
# the iters changes. the 9 iterations can be run all at once but is more computationally
# intensive. The data files below combine the results from iterations 1 to 4 and 5 to 9.
#########################################################################################
for(i in (1:n_iters)){
  cat("Repetition", i, '====================')
  eval(oneiteration)
}

######## Saving PLFAM Results ##########
date = format(Sys.time(), "%m-%d-%Y")
### Reps 1 to 4
# save(wpe2all, file = paste("OUT_PLFAM_MSPE_MIDW_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(wpe2all_mse_perc, file = paste("OUT_PLFAM_MSPE_PERC_MIDW_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(wpe2all_mse_var, file = paste("OUT_PLFAM_MSPE_VAR_MIDW_incl_irrig_reps1to4",  date, ".RData", sep = ""))
# save(wpe2all_mse_mean_abs, file = paste("OUT_PLFAM_MSPE_MeanAbs_MIDW_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(wpe2all_mse_mean_sqr, file = paste("OUT_PLFAM_MSPE_MeanSqr_MIDW_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(plfam_train_time, file = paste("OUT_PLFAM_TrainRunTime_incl_irrig_reps1to4", date, ".RData", sep = ""))

# Loading the data for reps 1 to 4 and setting new names for them so that later I can add the data for reps 5 to 9
# since they have the same name
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_PERC_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_VAR_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_MeanAbs_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_MeanSqr_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_TrainRunTime_incl_irrig_reps1to403-12-2024.RData")
wpe2all_reps1to4 = wpe2all
wpe2all_mse_perc_reps1to4 = wpe2all_mse_perc
wpe2all_mse_var_reps1to4 = wpe2all_mse_var
wpe2all_mse_mean_abs_reps1to4 = wpe2all_mse_mean_abs
wpe2all_mse_mean_sqr_reps1to4 = wpe2all_mse_mean_sqr
plfam_train_time_reps1to4 = plfam_train_time
### Reps 5 to 9
# save(wpe2all, file = paste("OUT_PLFAM_MSPE_MIDW_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(wpe2all_mse_perc, file = paste("OUT_PLFAM_MSPE_PERC_MIDW_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(wpe2all_mse_var, file = paste("OUT_PLFAM_MSPE_VAR_MIDW_incl_irrig_reps5to9",  date, ".RData", sep = ""))
# save(wpe2all_mse_mean_abs, file = paste("OUT_PLFAM_MSPE_MeanAbs_MIDW_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(wpe2all_mse_mean_sqr, file = paste("OUT_PLFAM_MSPE_MeanSqr_MIDW_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(plfam_train_time, file = paste("OUT_PLFAM_TrainRunTime_incl_irrig_reps5to9", date, ".RData", sep = ""))
# Loading the data for reps 5 to 9 
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_PERC_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_VAR_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_MeanAbs_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_MSPE_MeanSqr_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_PLFAM_TrainRunTime_incl_irrig_reps5to903-14-2024.RData")
wpe2all_reps5to9= wpe2all
wpe2all_mse_perc_reps5to9 = wpe2all_mse_perc
wpe2all_mse_var_reps5to9 = wpe2all_mse_var
wpe2all_mse_mean_abs_reps5to9 = wpe2all_mse_mean_abs
wpe2all_mse_mean_sqr_reps5to9 = wpe2all_mse_mean_sqr
plfam_train_time_reps5to9 = plfam_train_time

######## Saving FLM Results ##########
# save(wpe2linearall, file = paste("OUT_FLM_MSPE_MIDW_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(wpe2linearall_mse_perc, file = paste("OUT_FLM_MSPE_PERC_MIDW_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(wpe2linearall_mse_var, file = paste("OUT_FLM_MSPE_VAR_MIDW_incl_irrig_reps1to4",  date, ".RData", sep = ""))
# save(wpe2linearall_mse_mean_abs, file = paste("OUT_FLM_MSPE_MeanAbs_MIDW_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(wpe2linearall_mse_mean_sqr, file = paste("OUT_FLM_MSPE_MeanSqr_MIDW_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(flm_train_time, file = paste("OUT_FLM_TrainRunTime_incl_irrig_reps1to4", date, ".RData", sep = ""))
# save(totalweight2, file = paste("OUT_TotalWeight_incl_irrig_reps1to4", date, ".RData", sep = ""))

# Loading the data for reps 1 to 4 and setting new names for them so that later I can add the data for reps 5 to 9
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_PERC_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_VAR_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_MeanAbs_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_MeanSqr_MIDW_incl_irrig_reps1to403-12-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_TrainRunTime_incl_irrig_reps1to403-12-2024.RData")
wpe2linearall_reps1to4 = wpe2linearall
wpe2linearall_mse_perc_reps1to4 = wpe2linearall_mse_perc
wpe2linearall_mse_var_reps1to4 = wpe2linearall_mse_var
wpe2linearall_mse_mean_abs_reps1to4 = wpe2linearall_mse_mean_abs
wpe2linearall_mse_mean_sqr_reps1to4 = wpe2linearall_mse_mean_sqr
flm_train_time_reps1to4 = flm_train_time

# Loading the total weight for reps 1 to 4
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/OUT_TotalWeight_incl_irrig_reps1to403-12-2024.RData")
totalweight2_reps1to4 = totalweight2

######## Loading the data for reps 5 to 9
# save(wpe2linearall, file = paste("OUT_FLM_MSPE_MIDW_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(wpe2linearall_mse_perc, file = paste("OUT_FLM_MSPE_PERC_MIDW_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(wpe2linearall_mse_var, file = paste("OUT_FLM_MSPE_VAR_MIDW_incl_irrig_reps5to9",  date, ".RData", sep = ""))
# save(wpe2linearall_mse_mean_abs, file = paste("OUT_FLM_MSPE_MeanAbs_MIDW_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(wpe2linearall_mse_mean_sqr, file = paste("OUT_FLM_MSPE_MeanSqr_MIDW_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(flm_train_time, file = paste("OUT_FLM_TrainRunTime_incl_irrig_reps5to9", date, ".RData", sep = ""))
# save(totalweight2, file = paste("OUT_TotalWeight_incl_irrig_reps5to9", date, ".RData", sep = ""))
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_PERC_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_VAR_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_MeanAbs_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_MSPE_MeanSqr_MIDW_incl_irrig_reps5to903-14-2024.RData")
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/Midwest/OUT_FLM_TrainRunTime_incl_irrig_reps5to903-14-2024.RData")
wpe2linearall_reps5to9 = wpe2linearall
wpe2linearall_mse_perc_reps5to9 = wpe2linearall_mse_perc
wpe2linearall_mse_var_reps5to9 = wpe2linearall_mse_var
wpe2linearall_mse_mean_abs_reps5to9 = wpe2linearall_mse_mean_abs
wpe2linearall_mse_mean_sqr_reps5to9 = wpe2linearall_mse_mean_sqr
flm_train_time_reps5to9 = flm_train_time

# Loading the total weight for reps 1 to 4
load("Data_Extraction_cleaning/RdataFiles/PLFAM_FLM_Results/OUT_TotalWeight_incl_irrig_reps5to903-14-2024.RData")
totalweight2_reps5to9 = totalweight2

################# Combining Reps 1 to 4 with reps 5 to 9 for Midwest data ###################
######## PLFAM ##########
wpe2all_9iters = rbind(wpe2all_reps1to4, wpe2all_reps5to9[5:9,])
wpe2all_mse_perc_9iters = rbind(wpe2all_mse_perc_reps1to4, wpe2all_mse_perc_reps5to9[5:9,])
wpe2all_mse_var_9iters = rbind(wpe2all_mse_var_reps1to4, wpe2all_mse_var_reps5to9[5:9,])
wpe2all_mse_mean_abs_9iters = rbind(wpe2all_mse_mean_abs_reps1to4, wpe2all_mse_mean_abs_reps5to9[5:9,])
wpe2all_mse_mean_sqr_9iters = rbind(wpe2all_mse_mean_sqr_reps1to4, wpe2all_mse_mean_abs_reps5to9[5:9,])
plfam_train_time_9iters = rbind(plfam_train_time_reps1to4, plfam_train_time_reps5to9[5:9,])

######## FLM ##########
wpe2linearall_9iters = rbind(wpe2linearall_reps1to4, wpe2linearall_reps5to9[5:9,])
wpe2linearall_mse_perc_9iters = rbind(wpe2linearall_mse_perc_reps1to4, wpe2linearall_mse_perc_reps5to9[5:9,])
wpe2linearall_mse_var_9iters = rbind(wpe2linearall_mse_var_reps1to4, wpe2linearall_mse_var_reps5to9[5:9,])
wpe2linearall_mse_mean_abs_9iters = rbind(wpe2linearall_mse_mean_abs_reps1to4, wpe2linearall_mse_mean_abs_reps5to9[5:9,])
wpe2linearall_mse_mean_sqr_9iters = rbind(wpe2linearall_mse_mean_sqr_reps1to4, wpe2linearall_mse_mean_abs_reps5to9[5:9,])
flm_train_time_9iters = rbind(flm_train_time_reps1to4, flm_train_time_reps5to9)

# Total weight array 
totalweight2_9iters = rbind(totalweight2_reps1to4, totalweight2_reps5to9[5:9,])

############# PLFAM Prediction Errors ################
plfam_mse <- apply(wpe2all_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # WPE of PLFAM(joint)
plfam_mse_perc <- apply(wpe2all_mse_perc_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # PLFAM MSPE perc
plfam_mse_var <-  apply(wpe2all_mse_var_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # PLFAM MSPE var
plfam_mse_mean_abs <- apply(wpe2all_mse_mean_abs_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # PLFAM MSPE mean abs
plfam_mse_mean_sqr <- apply(wpe2all_mse_mean_sqr_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # PLFAM MSPE mean square
# Creating a dataframe to store all the PLFAM MSPE's
plfam_mspe_midw = data.frame(matrix(0, nr = 9, nc = 5))
colnames(plfam_mspe_midw) = c("test_MSE", "test_mse_perc", "test_mse_var", "test_mse_mean_abs", "test_mse_mean_sqr")
plfam_mspe_midw[,1] = plfam_mse
plfam_mspe_midw[,2] = plfam_mse_perc
plfam_mspe_midw[,3] = plfam_mse_var
plfam_mspe_midw[,4] = plfam_mse_mean_abs
plfam_mspe_midw[,5] = plfam_mse_mean_sqr
plfam_mspe_midw
date = format(Sys.time(), "%m-%d-%Y")
# Saving PLFAM Results
#save(plfam_mspe_midw, file = paste("PLFAM_MSPE_MIDW_incl_irrig_FINAL", date, ".RData", sep = ""))

############# FLM Prediction Errors ################
flm_mse <- apply(wpe2linearall_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # WPE of PLFAM(joint)
flm_mse_perc <- apply(wpe2linearall_mse_perc_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # PLFAM MSPE perc
flm_mse_var <-  apply(wpe2linearall_mse_var_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # PLFAM MSPE var
flm_mse_mean_abs <- apply(wpe2linearall_mse_mean_abs_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # PLFAM MSPE mean abs
flm_mse_mean_sqr <- apply(wpe2linearall_mse_mean_sqr_9iters, 1, sum)/apply(totalweight2_9iters, 1, sum) # PLFAM MSPE mean square
# Creating a dataframe to store all the PLFAM MSPE's
flm_mspe_midw = data.frame(matrix(0, nr = 9, nc = 5))
colnames(flm_mspe_midw) = c("test_MSE", "test_mse_perc", "test_mse_var", "test_mse_mean_abs", "test_mse_mean_sqr")
flm_mspe_midw[,1] = flm_mse
flm_mspe_midw[,2] = flm_mse_perc
flm_mspe_midw[,3] = flm_mse_var
flm_mspe_midw[,4] = flm_mse_mean_abs
flm_mspe_midw[,5] = flm_mse_mean_sqr
flm_mspe_midw
date = format(Sys.time(), "%m-%d-%Y")
# Saving PLFAM Results
#save(flm_mspe_midw, file = paste("FLM_MSPE_MIDW_incl_irrig_FINAL", date, ".RData", sep = ""))















