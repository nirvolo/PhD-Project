#' Making predictions from a fitted quantile spatially varying coefficient model
#' @param spline_coefs The estimated spline coefficients from the training model
#' \cr
#' @param Xpred The design matrix for prediction.
#' \cr
#' @param Spred The coordinates for prediction. If train this is Null
#' \cr
#' @param triang The triangulation of the spatial domain
#' \cr
#' @param B The basis function matrix
#' \cr
#' @param Q2 The Q2 matrix from the QR decomposition of the smoothness matrix H
#' \cr
#' @param train 
#' @return A vector of predicted response is returned.

predict_svqfm = function(spline_coefs, Xpred, Spred = NULL, triang = spat_tri, B = NULL, Q2 = NULL, train,
                         d = 3, r = 1){
  if(!is.matrix(Xpred)){
    warning("The explanatory variable, Xpred, should be a matrix.")
    Xpred = as.matrix(Xpred)
  }
  if(!is.null(Spred)){ # Spred is not required for training
    if(!is.matrix(Spred)){
      warning("The coordinates, Spred, should be a matrix.")
      Spred = as.matrix(Spred)
    }
  }

  if(train){ # Predictions for training data
    Bstar = B%*%Q2
    eta = Bstar%*%spline_coefs 
    ypred = rowSums(Xpred*eta)
    #cat("The length of ypred is", length(ypred), "\n")
  } else{ # predictions for testing data
    V = triang$V; Tr = triang$Tr
    # The basis functions for the training locations
    Basis.full = basis(V, Tr, d, r, Spred, Hmtx = FALSE, Kmtx = FALSE)
    ind.inside.pred = Basis.full$Ind.inside
    Bpred = Basis.full$B
    Bstar.pred = Bpred%*%Q2
    Xpred = as.matrix(Xpred[ind.inside.pred, ])
    # Calculating the estimated y values
    eta = Bstar.pred%*%spline_coefs
    ypred = rowSums(Xpred*eta)
  }
  return(list("yhat" = ypred, "eta" = eta))
}















