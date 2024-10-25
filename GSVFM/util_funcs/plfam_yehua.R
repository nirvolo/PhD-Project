########################
#### load libraries ####
########################
#install.packages("fda")
#install.packages("cosso")
#devtools::install_github("https://github.com/cran/cosso") another way to install cosso
library("fda")
library("cosso")


###################################
#### smoothing functional data ####
###################################
# for applying smoothing with automatic smoothing parameter selection
# using GCV

get.fd <- function(x, y, bs, loglam.lims=c(-15,5), lam.len=50, dfscale=1){
  # x           : a vector of time points (indices)
  # y           : a vector of observed values
  # bs          : basisfd object
  # loglam.lims : limits of the search of smoothing parameter (log-scale)
  # lam.len     : grid size of the search of smoothing parameter
  # dfscale     : scale parameter multiplied to the effective df of the model
  #               typically set as 1; for less wiggly model, set > 1

  loglam <- seq(loglam.lims[1], loglam.lims[2], len=lam.len)
  gcvs <- array(dim=lam.len)
  mse <- array(dim=lam.len)
  dfs <- array(dim=lam.len)
  n <- length(x)

  sobjs <- list()
  for (i in (1:lam.len)){
    lambdai <- exp(loglam[i])
    fP <- fdPar(bs, 2, lambdai)
    sobjs[[i]] <- smooth.basis(x, y, fP)
    gcvs[i] <- (sobjs[[i]]$SSE/n)/((max(0,1-sobjs[[i]]$df*dfscale/n))^2)
    mse[i] <- sobjs[[i]]$SSE/n
    dfs[i] <- sobjs[[i]]$df
  }
  i <- which.min(gcvs)
  if ((i==1)||(i==lam.len)){
    cat(c("warnings: selected smoothing parameters: i=",i,"\n"))
    warn <- T
  } else {
    warn <- F
  }
  return(list(sobj=sobjs[[i]], sobjs=sobjs, gcvs=gcvs, mses=mse, mse=mse[i], dfs=dfs,
              loglam=loglam, warn=warn))
}

# one smoothing parameter for all
# same time points
dat2fd2 <- function(argvals, y, rangeval, bs=NULL, loglam.lims=c(-15,5), lam.len=50,
                   dfscale=1, traceit=F, norder=4, nk=50)
{
  # argvals  :  see smooth.basis from fda library
  # y         : see smooth.basis from fda library
  # rangeeval : lower and upper end of time interval

  if (is.null(bs)){
    norder <- norder
    nk <- nk
    #x <- sort(argvals); n <- length(x); nk <- 20
    #times <- quantile(x[2:(n-1)], seq(0, 1, length=nk)) # knots
    times <- seq(rangeval[1], rangeval[2], length=nk)
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
  gcvId <- which.min(gcvs)
  sobj <- sobjs[[gcvId]]

  return(list(fd=sobj$fd, gcvs=gcvs, mmse=mmse[gcvId], mmses=mmse, dfs=dfs,
              lambda=exp(loglam[gcvId])))
}

# projecting new fd obj (represented by the same basis as the fdpca) to harmonics
# pcafd: pca object
# fd: fd object using the same basis
proj.fpca <- function(pcafd, fdobj){
  if (!(pcafd$harmonics$basis == fdobj$basis)){
    stop("basis mismatched.")
  }

  coef <- fdobj$coefs
  coefd <- dim(coef)
  nrep <- coefd[2]
  nharm <- dim(pcafd$harmonics$coefs)[2]

  # de-mean by the pcafd$meanfd
  fdobj$coefs <- fdobj$coefs - as.vector(pcafd$meanfd$coefs)

  harmscr <- array(0, c(nrep, nharm))
  coefarray <- fdobj$coefs
  harmcoefarray <- pcafd$harmonics$coefs
  basisobj <- fdobj$basis
  fdobjj <- fd(as.matrix(coefarray), basisobj)
  harmfdj <- fd(as.matrix(harmcoefarray), basisobj)
  harmscr <- inprod(fdobjj, harmfdj)
  return (harmscr)
}

proj.mfpca <- function(pcafd, fdobj){
  if (!(pcafd$harmonics$basis == fdobj$basis)){
    stop("basis mismatched.")
  }

  coef <- fdobj$coefs
  coefd <- dim(coef)
  nvar <- coefd[3]
  nrep <- coefd[2]
  nharm <- dim(pcafd$harmonics$coefs)[2]

  # de-mean by the pcafd$meanfd
  for (i in (1:nvar)){
    fdobj$coefs[,,i] <- fdobj$coefs[,,i] - pcafd$meanfd$coefs[,,i]
  }

  harmscr <- array(0, c(nrep, nharm, nvar))
  coefarray <- fdobj$coefs
  harmcoefarray <- pcafd$harmonics$coefs
  basisobj <- fdobj$basis
  for (j in 1:nvar){
    fdobjj <- fd(as.matrix(coefarray[, , j]), basisobj)
    harmfdj <- fd(as.matrix(harmcoefarray[, , j]), basisobj)
    harmscr[, , j] <- inprod(fdobjj, harmfdj)
  }
  return (harmscr)
}


########################
#### misc functions ####
########################

#### generate zeta
# fdobj: fd object
# nharm: number of harmonics (if NULL, thresh will be used to determine nharm). This is the number
# of principal components to compute.
# This function performs FPCA on the smoothed data X1 and X2 using the function
# pca.fd and then calculates zeta based on the principal component scores of X1 and X2.
get.zeta <- function(fdobj, nharm=NULL, thresh=0.999, centerfns=T){
  fpcaobj <- pca.fd(fdobj, nharm=fdobj$basis$nbasis, centerfns=centerfns)
  if (is.null(nharm))
    # ***It was fcumsum here instead of cumsum but I think it should be 
    # cumsum since I can't find a function like this. Also, in get.zeta2,
    # it's also cumsum***
    nharm <- sum(cumsum(fpcaobj$varprop) <= thresh)

  if (length(dim(fpcaobj$scores))<=2){
    score <- fpcaobj$scores
  } else {
    score <- apply(fpcaobj$scores, c(1,2), sum)
  }
  zeta <- pnorm(score[,1:nharm, drop=F], sd=rep(sqrt(fpcaobj$values[1:nharm]),
                                                each=nrow(score)))

  return(list(fpca=fpcaobj, nharm=nharm, zeta=zeta, centerfns=centerfns))
}


get.zeta2 <- function(nfobj, fdobj, nharm=NULL, thresh=0.999, centerfns=T){
  fpcaobj <- pca.fd(fdobj, nharm=fdobj$basis$nbasis, centerfns=centerfns)
  if (is.null(nharm))
    nharm <- sum(cumsum(fpcaobj$varprop) <= thresh)

  if (length(dim(fpcaobj$scores))<=2){
    score <- proj.fpca(fpcaobj, nfobj)
  } else {
    score <- proj.mfpca(fpcaobj, nfobj)
    score <- apply(score, c(1,2), sum)
  }
  zeta <- pnorm(score[,1:nharm, drop=F], sd=rep(sqrt(fpcaobj$values[1:nharm]),
                                                each=nrow(score)))
  return(list(fpca=fpcaobj, nharm=nharm, zeta=zeta))
}

# combine fd obj to a multivariate fd obj
# assuming same basis
combine.fd <- function(fdobj1, fdobj2){
  bscoef <- array(dim=c(dim(fdobj1$coefs),2))
  bscoef[,,1] <- fdobj1$coefs
  bscoef[,,2] <- fdobj2$coefs
  fdobj <- fd(bscoef, fdobj1$basis)
  return(fdobj)
}

############################################
#### fitting COSSO with parametric part ####
############################################
# from cosso packages: bigGram
# with weights

# from cosso package
solve.singular <- function (A, b) 
{
  solution = tryCatch(solve(A, b), error = function(x) NA)
  if (is.na(solution[1])) {
    solution = solve(A + 1e-07 * diag(nrow(A)), b)
  }
  return(solution)
}

# partial spline fitting (with weights)
psspline <- function(y, U1, R, phi, weight, kappa0=1)
{
  # y : response vector
  # U1 : design matrix of the parametric part (including the 1 vector)
  # R : 3d Gram
  # weight: weight vector: loss (1/n)\sum^n_{i=1} weight_i (y_i-f_i)^2 (req: weights sum to 1)
  # phi : the weights of the each component of R (along the third dimension)

  W <- diag(weight)
  s <- dim(R)[3]
  n <- dim(R)[1]
  p <- ncol(U1)
  RR <- matrix(0, nr=n, nc=n)
  for (j in (1:s)) RR <- RR + phi[j]*R[,,j]
  bigX <- cbind(U1, RR)  # big design matrix
  XX <- t(bigX) %*% W %*% bigX/n + kappa0*bdiag(matrix(0,nr=p,nc=p), RR)
  XXiX <- solve.singular(XX, t(bigX)) %*% W
  thetaa <- XXiX %*% y / n
  S <- bigX %*% XXiX / n         # smoothing matrix
  #return(as.vector(solve(XX, t(bigX)%*%y/n)))
  return(list(thetaa=thetaa, S=S, rss=sum(weight*(y-S%*%y)^2), weight=weight))
}

# find the best kappa for partial spline by minimizing the GCV
# ref. (3.34) of Gu13
psspline.gcv <- function(y, U1, R, phi, weight, kap=NULL)
{
  # argument definitions similar to psspline
  # kap : a vector of kappa candiates

  if (is.null(kap)) kap <- 2^seq(-10, -22, -0.75) # from cosso package

  # preparation
  W <- diag(weight)
  nkap <- length(kap)
  s <- dim(R)[3]
  n <- dim(R)[1]
  p <- ncol(U1)
  RR <- matrix(0, nr=n, nc=n)
  for (j in (1:s)) RR <- RR + phi[j]*R[,,j]
  bigX <- cbind(U1, RR)  # big design matrix
  bigXtXn <- t(bigX) %*% W %*% bigX / n

  # compute gcv
  gcv <- array(dim=nkap)
  for (k in (1:nkap)){
    XX <- bigXtXn + kap[k]*bdiag(matrix(0,nr=p,nc=p), RR)
    S <- bigX %*% solve.singular(XX, t(bigX)) %*% W / n         # smoothing matrix
    gcv[k] <- sum(weight*(y-S%*%y)^2)/(1-mean(diag(S)))^2
    # GCV (equivalent to (3.34) of Gu13)
  }

  jj <- which.min(gcv)
  if ((jj==1)||(jj==nkap)) warnings("psspline.gcv: need larger range of kappas:
                                    selected", jj, "\n")
  return(list(bestkap=kap[jj], gcv=gcv, kap=kap, weight=weight))
}


# two step fitting
twostepfit <- function(y, U1, R, kappa0, G, weight)
{
  # argument definitions similar to psspline
  # kappa0 : algorithm tuning parameter
  # G : tuning parameter of the the partial linear COSSO 
  # weight: weight vector: loss (1/n)\sum^n_{i=1} weight_i (y_i-f_i)^2 (req. weight sum to 1)
  
  p <- ncol(U1)
  n <- nrow(U1)
  s <- dim(R)[3]

  # partial spline
  thetaa <- psspline(y, U1, R, phi=rep(1,s), weight=weight, kappa0=kappa0)$thetaa
  theta <- thetaa[1:p]; a <- thetaa[-(1:p)]

  # QP
  D <- matrix(nr=n, nc=s)
  W <- diag(weight)
  for (j in (1:s)) D[,j] <- R[,,j] %*% a
  Dmat <- 2 * t(D) %*% W %*% D / n
  dvec <- 2 * t(D) %*% (W %*% (y - U1 %*% theta) - 0.5 * n * kappa0 * a) / n
  Amat <- t(rbind(diag(s), rep(-1, s)))
  bvec <- c(rep(0,s), -G)
  #phi <- solve.QP(Dmat, dvec, Amat, bvec)$solution
  K <- mean(abs(Dmat));
  Dmat <- Dmat/K; dvec <- dvec/K # to avoid scale issue
  
  ans <- try(solve.QP(Dmat, dvec, Amat, bvec), silent=T)
  add.to.diag <- sum(diag(Dmat)/ncol(Dmat))*1e-10
  while (is.character(ans)){
    cat(".")
    add.to.diag <- add.to.diag*10
    ans <- try(solve.QP(Dmat+diag(add.to.diag, nrow(Dmat)), dvec, Amat, bvec), silent=T)
  }
  phi <- ans$solution
  phi[phi<1e-08] <- 0

  # partial spline
  res <- psspline(y, U1, R, phi=phi, weight=weight, kappa0=kappa0)
  theta <-res$thetaa[1:p]; a <- res$thetaa[-(1:p)]

  return(list(theta=theta, a=a, phi=phi, S=res$S, rss=res$rss, weight=weight))
}

# fitting partial linear COSSO with BIC/GCV for choosing the tuning parameter
# note GCV is used for tuning the kappa0 (algorithm parameter)
plcosso <- function(y, U, zeta, weight=rep(1,length(y)), traceit=F, fine=T, nfine=20, sig2=NULL,
                    G.gcv=F, dfscale=1)
{
  # y : response
  # U : design matrix (vector of 1 will be appended by this function) (linear part)
  # zeta: design matrix (cosso part)
  # weight: weight vector: loss (1/n)\sum^n_{i=1} weight_i (y_i-f_i)^2
  # sig2: true variance of the observational noise of response; if not NULL, a
  #       different bic formula will be used
  # G.gcv: if it is false, G is chosen by BIC
  # dfscale: scale parameter of the degree of freedom. if dfscale>1, we have oversmoothing
  

  if (!(prod(weight>=0))) stop("some weights < 0")
  weight <- weight/sum(weight)*length(y)

  if (traceit) cat("preparing...\n")
  if (is.null(U)){
    if (traceit) cat("U is NULL: no parametric term\n")
    U1 <- matrix(1, nr=length(y), nc=1)
  } else {
    U1 <- cbind(1, U)
  }
  p <- ncol(U1)
  n <- nrow(U1)
  s <- ncol(zeta) # s
  R <- bigGram(zeta, zeta) # 3d Gram (n * n * s)
  if (traceit) cat("setting algorithm parameter kappa0 (using GCV)...\n")
  kappa0 <- psspline.gcv(y, U1, R, phi=rep(1,s), weight=weight)$bestkap   # look for kappa0
  G <- 0.2  # start at 0.2
  count <- 0; Grid <- NULL;
  theta <- NULL; a <- NULL; phi <- NULL; rss <- NULL; dfs <- NULL # not efficient to let the R objects grow
  tmpphi <- rep(0,s)
  if (traceit) cat("begin search of best G (using BIC/GCV)\n")
  while ((sum(tmpphi> 1e-07) < s) && (count <= ifelse(s<=15, floor(2*s), s))){
    # condition similar to cosso package
    if (traceit) cat("    G:", G)
    count <- count + 1
    Grid <- c(Grid, G)
    res <- twostepfit(y, U1, R, kappa0, G, weight)
    theta <- cbind(theta, res$theta); a <- cbind(a, res$a)
    phi <- cbind(phi, res$phi); dfs <- c(dfs, sum(diag(res$S))); rss <- c(rss, res$rss)
    tmpphi <- res$phi
    if (count < 10) G <- G + 0.25
    else if ((10 <= count) && (count < 16)) G <- G + 0.5
    else if ((16 <= count) && (count < 20)) G <- G + 1
    else G <- G + 2
    if (traceit) cat("\tnumber of nonzero components of phi:", sum(tmpphi>1e-07), "\n")
  }
  if (is.null(sig2)){
    bic <- n * log(rss/n) + dfscale * dfs * log(n)
    #bic <- n * log(rss/n) + apply(phi,2,function(x){sum(x>0)}) * log(n)
  } else {
    bic <- rss/sig2 + dfscale * dfs * log(n) 
  }
  gcv <- rss / (1- dfscale * dfs/n)^2
  if (G.gcv){
    whbestG <- which.min(gcv)
  } else {
    whbestG <- which.min(bic)
  }


  if (fine && (whbestG>1) && (whbestG<length(bic))){
    if (traceit) cat("begin fine-tuning\n")
    Grid2 <- seq(Grid[whbestG-1], Grid[whbestG+1], len=nfine)
    for (i in (1:nfine)){
      G <- Grid2[i]
      if (traceit) cat("    G:", G)
      Grid <- c(Grid, G)
      res <- twostepfit(y, U1, R, kappa0, G, weight)
      theta <- cbind(theta, res$theta); a <- cbind(a, res$a)
      phi <- cbind(phi, res$phi); dfs <- c(dfs, sum(diag(res$S))); rss <- c(rss, res$rss)
      if (traceit) cat("\tnumber of nonzero components of phi:", sum(res$phi>1e-07), "\n")
    }
    oo <- order(Grid); Grid <- Grid[oo]; theta <- theta[,oo, drop=F]; a <- a[, oo, drop=F]
    phi <- phi[,oo, drop=F]; dfs <- dfs[oo]; rss <- rss[oo];
    if (is.null(sig2)){
      bic <- n * log(rss/n) + dfscale * dfs * log(n)
      #bic <- n * log(rss/n) + apply(phi,2,function(x){sum(x>0)}) * log(n)
    } else {
      bic <- rss/sig2 + dfscale * dfs * log(n) 
    }
    gcv <- rss / (1- dfscale * dfs/n)^2
    if (G.gcv){
      whbestG <- which.min(gcv)
    } else {
      whbestG <- which.min(bic)
    }
  }

  if (traceit) cat("finishing...\n")
  btheta <- theta[,whbestG]; ba <- a[, whbestG]; bphi <- phi[, whbestG];
  ind <- (bphi>1e-07)
  plcossoobj <- list(btheta=btheta, ba=ba, bphi=bphi, ind=ind, whbestG=whbestG,
                     theta=theta, a=a, phi=phi, dfs=dfs, rss=rss, Grid=Grid, bic=bic, gcv=gcv,
                     U=U, zeta=zeta, weight=weight, kappa0=kappa0, G.gcv=G.gcv)
  class(plcossoobj) <- "plcosso"
  return(plcossoobj)
}

gen.groups <- function(n, nfold){
  leave.out <- trunc(n/nfold)
  o <- sample(1:n)
  groups <- vector("list", nfold)
  for (j in (1:(nfold-1))){
    jj <- (1+(j-1)*leave.out)
    groups[[j]] <- (o[jj:(jj+leave.out-1)])
  }
  groups[[nfold]] <- o[(1+(nfold-1)*leave.out):n]
  return(groups=groups)
}

# fitting partial linear COSSO with k-fold CV for choosing the tuning parameter
# note GCV is used for tuning the kappa0 (algorithm parameter)
plcosso.kcv <- function(y, U, zeta, weight=rep(1,length(y)), traceit=F, fine=T, nfine=20, nfold=5)
{
  # y : response
  # U : design matrix (vector of 1 will be appended by this function) (linear part)
  # zeta: design matrix (cosso part)
  # weight: weight vector: loss (1/n)\sum^n_{i=1} weight_i (y_i-f_i)^2
  

  if (!(prod(weight>=0))) stop("some weights < 0")
  weight <- weight/sum(weight)*length(y)

  if (traceit) cat("preparing...\n")
  if (is.null(U)){
    if (traceit) cat("U is NULL: no parametric term\n")
    U1 <- matrix(1, nr=length(y), nc=1)
  } else {
    U1 <- cbind(1, U)
  }
  p <- ncol(U1)
  n <- nrow(U1)
  s <- ncol(zeta) # s
  R <- bigGram(zeta, zeta) # 3d Gram (n * n * s)
  if (traceit) cat("setting algorithm parameter kappa0 (using GCV)...\n")
  kappa0 <- psspline.gcv(y, U1, R, phi=rep(1,s), weight=weight)$bestkap   # look for kappa0
  G <- 0.2  # start at 0.2
  count <- 0; Grid <- NULL;
  cvs <- NULL

  if (traceit) cat("begin search of best G (k-fold CV)\n")
  # form fold
  groups <- gen.groups(n, nfold)

  while ((count <= ifelse(s<=15, floor(2*s), s))){
    # condition similar to cosso package
    if (traceit) cat("    G:", G, "\n")
    count <- count + 1
    Grid <- c(Grid, G)
    cv <- 0
    for (j in (1:nfold)){
      ind <- groups[[j]]
      res <- twostepfit(y[-ind], U1[-ind,,drop=F], R[-ind,-ind,,drop=F], kappa0, G, weight[-ind])

      # prediction on the leave-out
      RR <- matrix(0, nr=length(ind), nc=n-length(ind))
      for (j in (1:s)) RR <- RR + res$phi[j] * R[ind,-ind,j]
      bigX <- cbind(U1[ind,,drop=F], RR)
      cv <- cv + sum(weight[ind]*(y[ind] - as.vector(bigX %*% c(res$theta, res$a)))^2)
    }
    cvs <- c(cvs, cv)
    if (count < 10) G <- G + 0.25
    else if ((10 <= count) && (count < 16)) G <- G + 0.5
    else if ((16 <= count) && (count < 20)) G <- G + 1
    else G <- G + 2
  }
  whbestG <- which.min(cvs)


  if (fine && (whbestG>1) && (whbestG<length(cvs))){
    if (traceit) cat("begin fine-tuning\n")
    Grid2 <- seq(Grid[whbestG-1], Grid[whbestG+1], len=nfine)
    for (i in (1:nfine)){
      G <- Grid2[i]
      if (traceit) cat("    G:", G, "\n")
      Grid <- c(Grid, G)
      cv <- 0
      for (j in (1:nfold)){
        ind <- groups[[j]]
        res <- twostepfit(y[-ind], U1[-ind,,drop=F], R[-ind,-ind,,drop=F], kappa0, G, weight[-ind])

        # prediction on the leave-out
        RR <- matrix(0, nr=length(ind), nc=n-length(ind))
        for (j in (1:s)) RR <- RR + res$phi[j] * R[ind,-ind,j]
        bigX <- cbind(U1[ind,,drop=F], RR)
        cv <- cv + sum(weight[ind]*(y[ind] - as.vector(bigX %*% c(res$theta, res$a)))^2)
      }
      cvs <- c(cvs, cv)
    }
    oo <- order(Grid); Grid <- Grid[oo]; cvs <- cvs[oo]
    whbestG <- which.min(cvs)
  }

  if (traceit) cat("finishing...\n")
  # Adding a print out to see what the optimal G is
  cat("The optimal G is", Grid[whbestG])
  res <- twostepfit(y, U1, R, kappa0, Grid[whbestG], weight)
  btheta <- res$theta; ba <- res$a; bphi <- res$phi;
  ind <- (bphi>1e-07)
  plcossoobj <- list(btheta=btheta, ba=ba, bphi=bphi, ind=ind, whbestG=whbestG,
                     Grid=Grid, cvs=cvs, U=U, zeta=zeta, weight=weight,
                     kappa0=kappa0)
  class(plcossoobj) <- "plcosso.kcv"
  return(plcossoobj)
}

# predict function for plcosso class
predict.plcosso <- function(obj, Unew, zetanew, whG=NULL)
{
  # whG: which G in the obj; if set as NULL, whG = obj$whbestG
  if (is.null(whG)) whG <- obj$whbestG
  n1 <- nrow(zetanew)
  n <- nrow(obj$zeta)
  s <- ncol(obj$zeta)
  #print(dim(zetanew))
  #print(dim(obj$zeta))
  R <- bigGram(zetanew, obj$zeta)
  RR <- matrix(0, nr=n1, nc=n)
  for (j in (1:s)) RR <- RR + obj$phi[j, whG] * R[,,j]
  bigX <- cbind(1, Unew, RR) # append vector of 1
  # cat("Dim of bigX", dim(bigX), "\n")
  # cat("Dim of whG", length(c(obj$theta[,whG], obj$a[,whG])), "\n")
  return(as.vector(bigX %*% c(obj$theta[,whG], obj$a[,whG])))
}

# predict function for plcosso class
predict.plcosso.kcv <- function(obj, Unew, zetanew)
{
  # whG: which G in the obj; if set as NULL, whG = obj$whbestG
  n1 <- nrow(zetanew)
  n <- nrow(obj$zeta)
  s <- ncol(obj$zeta)
  R <- bigGram(zetanew, obj$zeta)
  RR <- matrix(0, nr=n1, nc=n)
  for (j in (1:s)) RR <- RR + obj$bphi[j] * R[,,j]
  bigX <- cbind(1, Unew, RR) # append vector of 1
  return(as.vector(bigX %*% c(obj$btheta, obj$ba)))
}


# internal function for plcosso.valid
predict.valid <- function(Unew, zetanew, zeta, phi, theta, a)
{
  n1 <- nrow(zetanew)
  n <- nrow(zeta)
  s <- ncol(zeta)
  R <- bigGram(zetanew, zeta)
  RR <- matrix(0, nr=n1, nc=n)
  for (j in (1:s)) RR <- RR + phi[j] * R[,,j]
  bigX <- cbind(Unew, RR) # do not append vector of 1
  return(as.vector(bigX %*% c(theta, a)))
}

# for plotting component function, assuming no interaction terms
comps <- function(obj, ngrid=100,  whG=NULL)
{
  if (is.null(whG)) whG <- obj$whbestG
  n <- nrow(obj$zeta)
  s <- ncol(obj$zeta)
  z <- seq(0, 1, len=ngrid)
  zetanew <- matrix(z, nr=ngrid, nc=s)
  R <- bigGram(zetanew, obj$zeta)
  comfit <- matrix(0,nr=ngrid, nc=s)
  for (j in (1:s)) if (obj$ind[j]) comfit[,j] <- obj$phi[j, whG] * R[,,j] %*% obj$a[,whG]
  return(list(z=z, comfit=comfit))
}

comps.kcv <- function(obj, ngrid=100,  whG=NULL)
{
  n <- nrow(obj$zeta)
  s <- ncol(obj$zeta)
  z <- seq(0, 1, len=ngrid)
  zetanew <- matrix(z, nr=ngrid, nc=s)
  R <- bigGram(zetanew, obj$zeta)
  comfit <- matrix(0,nr=ngrid, nc=s)
  for (j in (1:s)) if (obj$ind[j]) comfit[,j] <- obj$bphi[j] * R[,,j] %*% obj$ba
  return(list(z=z, comfit=comfit))
}

# fitting partial linear COSSO using validation set for choosing the tuning parameter
# note GCV is used for tuning the kappa0 (algorithm parameter)
plcosso.valid <- function(y, U, zeta, weight,
                          yv, Uv, zetav, weightv,
                          traceit=F, fine=T, nfine=20, sig2=NULL,
                          G.gcv=F, dfscale=1)
{
  # y : response
  # U : design matrix (vector of 1 will be appended by this function) (linear part)
  # zeta: design matrix (cosso part)
  # weight: weight vector: loss (1/n)\sum^n_{i=1} weight_i (y_i-f_i)^2
  # sig2: true variance of the observational noise of response; if not NULL, a
  #       different bic formula will be used
  # G.gcv: if it is false, G is chosen by BIC
  # dfscale: scale parameter of the degree of freedom. if dfscale>1, we have oversmoothing
  

  if (!(prod(weight>=0))) stop("some weights < 0")
  weight <- weight/sum(weight)*length(y)
  weightv <- weightv/sum(weightv)*length(yv)

  if (traceit) cat("preparing...\n")
  if (is.null(U)){
    if (traceit) cat("U is NULL: no parametric term\n")
    U1 <- matrix(1, nr=length(y), nc=1)
    Uv <- matrix(1, nr=length(yv), nc=1)
  } else {
    U1 <- cbind(1, U)
    Uv <- cbind(1, Uv)
  }
  p <- ncol(U1)
  n <- nrow(U1)
  s <- ncol(zeta) # s
  R <- bigGram(zeta, zeta) # 3d Gram (n * n * s)
  if (traceit) cat("setting algorithm parameter kappa0 (using GCV)...\n")
  kappa0 <- psspline.gcv(y, U1, R, phi=rep(1,s), weight=weight)$bestkap   # look for kappa0
  G <- 0.2  # start at 0.2
  count <- 0; Grid <- NULL;
  theta <- NULL; a <- NULL; phi <- NULL; rss <- NULL; dfs <- NULL # not efficient to let the R objects grow
  wpe <- NULL;
  tmpphi <- rep(0,s)
  if (traceit) cat("begin search of best G (using BIC/GCV)\n")
  while ((sum(tmpphi> 1e-07) < s) && (count <= ifelse(s<=15, floor(2*s), s))){
    # condition similar to cosso package
    if (traceit) cat("    G:", G)
    count <- count + 1
    Grid <- c(Grid, G)
    res <- twostepfit(y, U1, R, kappa0, G, weight)
    pp <- predict.valid(Uv, zetav, zeta, res$phi, res$theta, res$a)
    wpe <- c(wpe, mean((yv-pp)^2*weightv))
    theta <- cbind(theta, res$theta); a <- cbind(a, res$a)
    phi <- cbind(phi, res$phi); dfs <- c(dfs, sum(diag(res$S))); rss <- c(rss, res$rss)
    tmpphi <- res$phi
    if (count < 10) G <- G + 0.25
    else if ((10 <= count) && (count < 16)) G <- G + 0.5
    else if ((16 <= count) && (count < 20)) G <- G + 1
    else G <- G + 2
    if (traceit) cat("\tnumber of nonzero components of phi:", sum(tmpphi>1e-07), "\n")
  }
  if (is.null(sig2)){
    bic <- n * log(rss/n) + dfscale * dfs * log(n)
    #bic <- n * log(rss/n) + apply(phi,2,function(x){sum(x>0)}) * log(n)
  } else {
    bic <- rss/sig2 + dfscale * dfs * log(n) 
  }
  gcv <- rss / (1- dfscale * dfs/n)^2
  whbestG <- which.min(wpe)


  if (fine && (whbestG>1) && (whbestG<length(bic))){
    if (traceit) cat("begin fine-tuning\n")
    Grid2 <- seq(Grid[whbestG-1], Grid[whbestG+1], len=nfine)
    for (i in (1:nfine)){
      G <- Grid2[i]
      if (traceit) cat("    G:", G)
      Grid <- c(Grid, G)
      res <- twostepfit(y, U1, R, kappa0, G, weight)
      pp <- predict.valid(Uv, zetav, zeta, res$phi, res$theta, res$a)
      wpe <- c(wpe, mean((yv-pp)^2*weightv))
      theta <- cbind(theta, res$theta); a <- cbind(a, res$a)
      phi <- cbind(phi, res$phi); dfs <- c(dfs, sum(diag(res$S))); rss <- c(rss, res$rss)
      if (traceit) cat("\tnumber of nonzero components of phi:", sum(res$phi>1e-07), "\n")
    }
    oo <- order(Grid); Grid <- Grid[oo]; theta <- theta[,oo, drop=F]; a <- a[, oo, drop=F]
    phi <- phi[,oo, drop=F]; dfs <- dfs[oo]; rss <- rss[oo];
    wpe <- wpe[oo]
    if (is.null(sig2)){
      bic <- n * log(rss/n) + dfscale * dfs * log(n)
      #bic <- n * log(rss/n) + apply(phi,2,function(x){sum(x>0)}) * log(n)
    } else {
      bic <- rss/sig2 + dfscale * dfs * log(n) 
    }
    gcv <- rss / (1- dfscale * dfs/n)^2
    whbestG <- which.min(wpe)
  }

  if (traceit) cat("finishing...\n")
  btheta <- theta[,whbestG]; ba <- a[, whbestG]; bphi <- phi[, whbestG];
  ind <- (bphi>1e-07)
  plcossoobj <- list(btheta=btheta, ba=ba, bphi=bphi, ind=ind, whbestG=whbestG,
                     theta=theta, a=a, phi=phi, dfs=dfs, rss=rss, Grid=Grid, bic=bic, gcv=gcv,
                     U=U, zeta=zeta, weight=weight, kappa0=kappa0, G.gcv=G.gcv, wpe=wpe)
  class(plcossoobj) <- "plcosso.valid"
  return(plcossoobj)
}



###########################################
#### ridge regression (for comparison) ####
###########################################
# find the best kappa for partial spline by minimizing GCV
# ref. (3.34) of Gu13
# only penalize coeff of R
myridge <- function(y, U, R, weight = rep(1,length(y)), kap=NULL)
{
  # argument definitions similar to psspline
  # kap : a vector of kappa candiates (regularization parameter)

  if (is.null(U)) U1 <- matrix(1, nr=length(y), nc=1)
  else U1 <- cbind(1,U)
  if (is.null(kap)) kap <- 2^seq(-10, -22, -0.75) # from cosso package


  # preparation
  W <- diag(weight)
  nkap <- length(kap)
  n <- length(y)
  p <- ncol(U1)
  s <- ncol(R)
  bigX <- cbind(U1, R)  # big design matrix
  bigXtXn <- t(bigX) %*% W %*% bigX / n

  # compute gcv
  gcv <- array(dim=nkap)
  for (k in (1:nkap)){
    XX <- bigXtXn + kap[k]*diag(c(rep(0,p), rep(1, s)))
    S <- bigX %*% solve.singular(XX, t(bigX)) %*% W / n         # smoothing matrix
    gcv[k] <- sum(weight*(y-S%*%y)^2)/(1-mean(diag(S)))^2
    # GCV (equivalent to (3.34) of Gu13)
  }

  jj <- which.min(gcv)
  if ((jj==1)||(jj==nkap)) warnings("need larger range of kappas: selected", jj, "\n")

  XX <- bigXtXn + kap[jj]*diag(c(rep(0,p), rep(1, s)))
  beta <- solve.singular(XX, t(bigX)) %*% W %*% y / n         # smoothing matrix
  XX <- bigXtXn
  beta0 <- solve.singular(XX, t(bigX)) %*% W %*% y / n         # smoothing matrix
  return(list(beta=as.vector(beta), beta0=as.vector(beta0), bestkap=kap[jj], gcv=gcv, kap=kap, weight=weight))
}
