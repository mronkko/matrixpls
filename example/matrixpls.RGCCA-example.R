# Runs the example from the rgcca function of the RGCCA package

if(require(RGCCA) & require(Matrix)) {
  
  data(Russett)
  X_agric =as.matrix(Russett[,c("gini","farm","rent")])
  X_ind = as.matrix(Russett[,c("gnpr","labo")])
  X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
  A = list(X_agric, X_ind, X_polit)
  #Define the design matrix (output = C)
  C = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
  
  tau <- 0
  
  result.rgcca = rgcca(A, C, tau =rep(tau,3), scheme = "factorial", scale = TRUE)
  
  W.rgcca <- as.matrix(bdiag(lapply(result.rgcca$a, function(x){matrix(x,nrow=1)})))
  
  # Do the same with matrixpls
  
  W.mod <- (W.rgcca != 0) *1
  
  S <- cov(do.call(cbind,A))
  W.matrixpls <- weight.pls(S, list(inner = C,
                                    reflective = t(W.mod),
                                    formative = matrix(0,nrow(W.mod), ncol(W.mod))),
                            W.mod = W.mod,
                            innerEstimator = inner.factor,
                            outerEstimators = outer.RGCCA, tau = tau)
  
  print(W.rgcca)
  print(W.matrixpls)
  
  # Try a different tau value
  
  tau2 <- .5
  result.rgcca2 <- rgcca(A, C, tau =rep(tau2,3), scheme = "factorial", scale = TRUE)
  
  W.matrixpls2 <- weight.pls(S, list(inner = C,
                                    reflective = t(W.mod),
                                    formative = matrix(0,nrow(W.mod), ncol(W.mod))),
                            W.mod = W.mod,
                            innerEstimator = inner.factor,
                            outerEstimators = outer.RGCCA, tau = tau2)
  
  # Do we get perfect correlations for the composites
  cor(cbind(do.call(cbind,result.rgcca$Y), 
            do.call(cbind,A) %*% t(W.matrixpls),
      do.call(cbind,result.rgcca2$Y), 
      do.call(cbind,A) %*% t(W.matrixpls2))
  )
}
