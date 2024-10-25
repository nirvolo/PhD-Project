#############################
#### libraries and codes ####
#############################
library(MASS)
library(dplyr)
library(tidyr)
source("GSVFM/util_funcs/plfam_yehua.R")
source("GSVFM/PLFAM_FLM_ExtractID_Kansas.R")


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
    temp <-spread(fundat[, c("Year", "DateI", "Id", "TMAX")], DateI, TMAX)
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
    temp <-spread(fundat[, c("Year", "DateI", "Id", "TMAX")], DateI, TMAX)
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

create.year.dummy <- function(y, uyears){
  p <- length(uyears)
  out <- matrix(0,nr=length(y), nc=p-1)
  for (i in (1:(p-1))){
    out[y==uyears[i],i] <- 1
  }
  return(out)
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
###### Loading the data #######
load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcRegdat.RData")
load("Data_Extraction_cleaning/RdataFiles/NewKansasFromMidw_PreProcFundat.RData")
regdat0 = kns_nonfd_new
fundat0 = kns_fd_new

n_iters = 9 # number of repetitions
n_folds <- 5 # number of folds in validations
DE_MEAN_RESP = T
llam.lims = c(-15.5,-8)


# define expression for tryCatch
oneiteration <-expression({
  set.seed(12876+i*123)
  for(vv in 1:n_folds){
    cat("\n","Fold", vv, "\n")
    cat('------------', "\n")
    # Each column of test_inds represents one of the 45 folds
    Id.valid <- na.omit(test_inds[,(i-1)*n_folds + vv])
    ####### Splitting the data into training and testing sets #########
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
    for (yr in unique(regdat2$Year)) {
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
    U1 <- cbind(regdat1$irrig_prop, regdat1$avgPRCP,
                regdat1$irrig_prop*regdat1$avgPRCP)
    weight1 <- regdat1$Area
    
    #### joint #####
    # check id.TMIN and id.TMAX
    if (!any(id.TMIN!=id.TMAX)){
      # both sobj.TMAX and sobj.TMIN are constructed on the same basis
      fd.joint <- combine.fd(sobj.TMAX$fd, sobj.TMIN$fd) # because id.TMIN==id.TMAX
      zobj.joint <- get.zeta(fd.joint)
      xi.joint1 <- zobj.joint$zeta[getId(regdat1$Id, id.TMAX),] # because id.TMIN==id.TMAX
      fit2 <- plcosso.kcv(y1, U1, xi.joint1, weight=weight1/sum(weight1), traceit=T)
      fit2.linear <- myridge(y1, U1, xi.joint1, weight=weight1/sum(weight1))
    }

    ##########################
    #### prediction error ####
    ##########################
    if(!DE_MEAN_RESP) y2 <- regdat2$Yield
    else y2 = regdat2$de_meaned_Yield
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
    pred2 <- predict.plcosso.kcv(fit2, U2, xi.joint2)
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
  
    totalweight2[i, vv] <- sum(weight2)
}})

############################
#### run the experiment ####
############################
for(i in (1:n_iters)){
  cat("Repetition", i, '====================')
  eval(oneiteration)
}

#####################
#### get summary ####
#####################
############# PLFAM Prediction Errors ################
owpe2all <- apply(wpe2all, 1, sum)/apply(totalweight2, 1, sum) # WPE of PLFAM(joint)
owpe2all_mse_perc <- apply(wpe2all_mse_perc, 1, sum)/apply(totalweight2, 1, sum) # PLFAM MSPE perc
owpe2all_mse_var <-  apply(wpe2all_mse_var, 1, sum)/apply(totalweight2, 1, sum) # PLFAM MSPE var
owpe2all_mse_mean_abs <- apply(wpe2all_mse_mean_abs, 1, sum)/apply(totalweight2, 1, sum) # PLFAM MSPE mean abs
owpe2all_mse_mean_sqr <- apply(wpe2all_mse_mean_sqr, 1, sum)/apply(totalweight2, 1, sum) # PLFAM MSPE mean square
plfam_mspe = data.frame(matrix(0, nr = 9, nc = 5))
colnames(plfam_mspe) = c("test_MSE", "test_mse_perc", "test_mse_var", "test_mse_mean_abs", "test_mse_mean_sqr")
plfam_mspe[,1] <- owpe2all
plfam_mspe[,2] <- owpe2all_mse_perc
plfam_mspe[,3] <- owpe2all_mse_var
plfam_mspe[,4] <- owpe2all_mse_mean_abs
plfam_mspe[,5] <- owpe2all_mse_mean_sqr
#save(plfam_mspe, file = "NewKnsData_PLFAM_MSPE_withKFoldCV_03_11_2024.RData")

############ FLM Prediction Errors ##################
owpe2linearall <- apply(wpe2linearall, 1, sum)/apply(totalweight2, 1, sum) # WPE of FLM-Cov(joint)
owpe2linearall_mse_perc <- apply(wpe2linearall_mse_perc, 1, sum)/apply(totalweight2, 1, sum) # FLM MSPE percent
owpe2linearall_mse_var <- apply(wpe2linearall_mse_var, 1, sum)/apply(totalweight2, 1, sum) # FLM MSPE var
owpe2linearall_mse_mean_abs <- apply(wpe2linearall_mse_mean_abs, 1, sum)/apply(totalweight2, 1, sum) # FLM MSPE mean abs
owpe2linearall_mse_mean_sqr <- apply(wpe2linearall_mse_mean_sqr, 1, sum)/apply(totalweight2, 1, sum) # FLM MSPE mean square
# Creating a dataframe to store all the FLM MSPE's
flm_mspe = data.frame(matrix(0, nr = 9, nc = 5))
colnames(flm_mspe) = c("test_MSE", "test_mse_perc", "test_mse_var", "test_mse_mean_abs", "test_mse_mean_sqr")
flm_mspe[,1] = owpe2linearall
flm_mspe[,2] = owpe2linearall_mse_perc
flm_mspe[,3] = owpe2linearall_mse_var
flm_mspe[,4] = owpe2linearall_mse_mean_abs
flm_mspe[,5] = owpe2linearall_mse_mean_sqr
flm_mspe
#save(flm_mspe, file = "NewKnsData_FLM_MSPE_withKFoldCV_03_11_2024.RData")
