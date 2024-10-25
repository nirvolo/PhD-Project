#' Making predictions from a fitted generalized spatially varying coefficient model
#'
#' This function is used to make predictions of the generalized spatially varying coefficient models.
#'
#' @importFrom BPST basis
#'
#' @param mfit Fitted ``gsvcm" object.
#' \cr
#' @param Xpred The design matrix for prediction.
#' \cr
#' @param Spred The cooridinates for prediction.
#' \cr
#' @return A vector of predicted response is returned.


predict.gsvcm = function(mfit, Xpred, Spred){
  #print(dim(Xpred))
  #print(dim(Spred))
  if(!is.matrix(Xpred)){
    warning("The explanatory variable, Xpred, should be a matrix.")
    Xpred = as.matrix(Xpred)
  }
  if(!is.matrix(Spred)){
    warning("The coordinates, Spred, should be a matrix.")
    Spred = as.matrix(Spred)
  }

  family = mfit$family; linkinv = family$linkinv;
  
  
  #print(Spred)
  #print(unique(mfit$S))
  #print(identical(Spred, unique(mfit$S)))
  #print(dim(Spred))
  #print(Spred)
  #print(dim(mfit$S))
  #print(identical(Spred, mfit$S))
  if(identical(Spred, mfit$S)){
  # I changed the identical condition here because it gave me an error when I 
  # didn't use the unique number of locations in mfit$S. The reason it gave an error
  # is because in Spred there are only 100 observations but in mfit$S there are all 
  # the training observations (900 in 10-fold CV). This is because mfit$S is coming from 
  # the training data where there are repeated spatial locations in the location matrix.
  #if(identical(Spred, unique(mfit$S))){
      W = as.matrix(kr(mfit$X, mfit$B %*% mfit$Q2, byrow = TRUE))
      #print(dim(W))
      eta = W %*% as.vector(mfit$theta_hat)
      ypred = linkinv(eta)
    } else {
      V = mfit$V; Tr = mfit$Tr; d = mfit$d; r = mfit$r; Q2 = mfit$Q2
      
      Basis.full = basis(V, Tr, d, r, Spred, FALSE, FALSE)
      ind.inside.pred = Basis.full$Ind.inside
      Bpred = Basis.full$B
      #print(c("The dimension of Bpred is",dim(Bpred)))
      #print(c("The dimension of Q2 is",dim(Q2)))
      Xpred = as.matrix(Xpred[ind.inside.pred, ])
      #print(c("The dimension of Xpred is", dim(Xpred)))
    
      W = as.matrix(kr(Xpred, Bpred %*% Q2, byrow = TRUE)) 
      #print(W)
      #print(c("The dimension of the W matrix is",dim(W)))
      #print(c("The dimension of the theta_hat matrix is", dim(mfit$theta_hat)))
      #print(length(as.vector(mfit$theta_hat)))
      eta = W %*% as.vector(mfit$theta_hat)
      #print(c("This is eta", eta))
      ypred = linkinv(eta)
      #print(c("This is ypred",ypred))
    }
  return(ypred)
}
