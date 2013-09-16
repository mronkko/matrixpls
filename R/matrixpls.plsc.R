#
# This R code is contributed to matrixpls by Huang. It is a part of her dissertation.
##
# Huang, W. (2013). PLSe: Efficient Estimators and Tests for Partial Least Squares
# (Doctoral dissertation). University of California, Los Angeles.
#

# =========== matrixpls specific changes ===========



# =========== Original code by Huang ===========

#This R code is based on Dijkstra's PLSc_UCLA.m with addition of situations when 
#there is only one observed variable loading on a single phantom latent variable,
#so that this one observed variable can serve as both predictor or outcome later on
#when we estimate the structural model regression coefficients.

#define an index function 
PLSc.idx <- function(p,i) {
  (1+p[i]):p[i+1]
}

# Determination of the sign-factors, procedure not described in the paper
# Wold's adjacency concept is not used: all latent variables are taken to be adjacent.
#ab is number of blocks/factors, 
PLSc.sign <- function(ab,p,W,S) {
  s <- matrix(0,ab,ab) #to set up a factor corr matrix
  for (i in 1:(ab-1)) {
    for (j in (i+1):ab) {
      idxi <- PLSc.idx(p,i)
      idxj <- PLSc.idx(p,j)
      s[i,j] <- sign(t(W[idxi])%*%S[idxi,idxj]%*%W[idxj])
    }
  }
  s <- t(s)+s
  return(s)
}

# Normalization of the weight vectors, equation 6 of Dijkstra, April 7, 2011.
PLSc.normalize <- function(ab,p,W,S) {
  for (i in 1:ab) {
    idx <- PLSc.idx(p,i)
    if (length(idx) > 1) { # only for latent factors
      W[idx] <- W[idx]/sqrt(t(W[idx])%*%S[idx,idx]%*%W[idx])
    }
  }
  return(W)
}

# Run through equation (4) of Dijkstra, April 7, 2011.
PLSc.weights <- function(ab,p,Wold,s,S) {
  Wnew <- Wold
  for (i in 1:ab) {
    # Example: idxi=1:3 when i=1, p=c(0 3 6 7), note p is not the p user defined
    idxi <- PLSc.idx(p,i)
    if (length(idxi) > 1) { # if more than one indicator per latent factor, then update the weight.
    # when there is only one indicator per factor, length(idxi)=1, don't update, weight is 1.
      Wnew[idxi] <- 0
      for (j in 1:ab) {
        idxj <- PLSc.idx(p,j)
        Wnew[idxi] <- Wnew[idxi]+s[i,j]*(S[idxi,idxj]%*%matrix(Wold[idxj]))
      }
    }
  }
  return(Wnew)
}

PLSc <- function(S,p,tol=1e-6,tmax=500) {
  ab <- length(p) #number of blocks
  ai <- sum(p) #total number of indicators
  p <- c(0,cumsum(p)) # column indices, e.g if p=c(3,3,1) then p=c(0,3,6,7)
  
  # Normalization of the first weight vectors
  Wold <- PLSc.normalize(ab,p,rep(1,ai),S) #The starting values for W are vectors proportional to(1,1,...,1)',
  # iteration to convergence or max iteration exceeded
  for (t in 1:tmax) {
    s <- PLSc.sign(ab,p,Wold,S) # Determination of the sign-factors
    Wnew <- PLSc.weights(ab,p,Wold,s,S) # weights
    Wnew <- PLSc.normalize(ab,p,Wnew,S) # normalize again
    Wnew <- abs(Wnew) # make loadings positive: see note 2 in Matlab code
    if (max(abs(Wnew-Wold)) < tol) { # test whether convergence has occurred
      break # if yes, break out of the iterations
    } else { # if not, keep on iterating
      Wold <- Wnew
    }
  }

  # End of iterations, determination of correlations, loadings and the quality of the proxies
  wdak <- Wnew # wdak contains the PLS mode A estimates for the weights.


  # Calculation of the correlations between the PLS mode A proxies, Rdak:
  Rdak <- matrix(0,ab,ab)
  for (i in 1:(ab-1)) {
    for (j in (i+1):ab) {
      idxi <- PLSc.idx(p,i)
      idxj <- PLSc.idx(p,j)
      Rdak[i,j] <- t(wdak[idxi])%*%S[idxi,idxj]%*%wdak[idxj]
    }
  }
  Rdak <- Rdak+t(Rdak)+diag(ab)
###may use Rdak as PLS estimates to compare with PLSc estimates R###


  # Determination of the correction factors, based on (11) of Dijkstra, April 7, 2011.
  c2 <- rep(1,ab)
  for (i in 1:ab) {
    idx <- PLSc.idx(p,i)
    if (length(idx) > 1) { # only for latent factors, no need to correct for the single indicator for the phantom LV
      c2[i] <- t(wdak[idx])%*%(S[idx,idx]-diag(diag(S[idx,idx])))%*%wdak[idx]
      c2[i] <- c2[i]/(t(wdak[idx])%*%(wdak[idx]%*%t(wdak[idx])-diag(diag(wdak[idx]%*%t(wdak[idx]))))%*%wdak[idx])
    }
  }
  c <- sqrt(c2)

  # Determination of consistent estimates of the loadings, see (13) of Dijkstra, April 7, 2011.
  L <- wdak
  for (i in 1:ab) {
    idx <- PLSc.idx(p,i)
    L[idx] <- c[i]*wdak[idx]
  }

  # Determination of the quality of the proxies, see (15) of Dijkstra, April 7, 2011.
  Q <- c2
  for (i in 1:ab) {
    idx <- PLSc.idx(p,i)
    Q[i] <- c2[i]*(t(wdak[idx])%*%wdak[idx])^2
  }

  # Determination of consistent estimates for the correlation between the
  # latent variables, see (15) and (16) of Dijkstra, April 7, 2011.
  R <- matrix(0,ab,ab)
  # R start with a matrix with all 0, then define the upper triangle of the corr matrix
  for (i in 1:(ab-1)) {
    for (j in (i+1):ab) {
      R[i,j] <- Rdak[i,j]/(sqrt(Q[i]*Q[j]))
    }
  }
  R <- R+t(R)+diag(ab) # to form the complete corr matrix of latent variables
  return(list(R0=Rdak,R=R,L=L,Q=Q,t=t))
} # R0 is PLS estimates for LV corr without Dijkstra's correction.



##Use Two Stage Least Squares to find the structural relationship between the laten factors
TSLS_findindex <- function(v) {
  I <- NULL
  for (j in 1:length(v)) {
    if (v[j] == 1) {
      I <- c(I,j)
    }
  }
  return(I)
}

TSLS_general <- function(R,IB,IC) {
  n <- nrow(IB) # number of etas
  x <- ncol(IC) # number of ksis
  B <- matrix(0,n,n)
  C <- matrix(0,n,x)
  H <- solve(R[1:x,1:x])
  
  for (i in 1:n) {
    I <- TSLS_findindex(IB[i,])
    J <- TSLS_findindex(IC[i,])
    HH <- rbind(cbind(R[x+I,1:x]%*%H%*%R[1:x,x+I],matrix(R[x+I,J],length(x+I),length(J))),
                cbind(matrix(R[J,x+I],length(J),length(x+I)),matrix(R[J,J],length(J),length(J))))
    HH <- ginv(HH)
    h <- rbind(R[x+I,1:x]%*%H%*%R[1:x,i+x],matrix(R[J,i+x]))
    theta <- HH%*%h
    B[i,I] <- theta[1:sum(IB[i,])]
    C[i,J] <- theta[(sum(IB[i,])+1):length(theta)]
  } 
  cov_zeta <- (diag(1,n)-B)%*%R[(x+1):(x+n),(x+1):(x+n)]%*%t(diag(1,n)-B)-C%*%R[(1:x),(1:x)]%*%t(C)
  return(list(B=B,C=C,cov=cov_zeta))
}



# takes covariance matrix return normal theory information matrix
FisherInfo <- function(V) {
  p <- nrow(V)
  q <- p*(p+1)/2
  invV <- solve(V) #solve returns inverse of V
  I <- matrix(0,q,q)

  row <- 0
  for (i in 1:p) {
    for (j in 1:i) {
      row <- row + 1
      col <- 0
      # derivative of V with respect to element (i,j)
      dV1 <- matrix(0,p,p)
      dV1[i,j] <- 1
      dV1[j,i] <- 1

      for (k in 1:p) {
        for (l in 1:k) {
          col <- col + 1
          # derivative of V wrt element (k,l)
          dV2 <- matrix(0,p,p)
          dV2[k,l] <- 1
          dV2[l,k] <- 1
          I[row,col] <- 0.5*sum(diag(invV%*%dV1%*%invV%*%dV2))
#          print(c(row,col))
        }
      }
    }
  }
  return(I)
}

# utilities

vech <- function(V) {
  d <- nrow(V)
  v <- rep(0,d*(d+1)/2)
  count <- 1
  for (i in 1:d) {
    for (j in 1:i) {
      v[count] <- V[i,j]
      count <- count+1
    }
  }
  return(v)
}


vech2 <- function(V) {
  d <- nrow(V)
  v <- rep(0,d*(d+1)/2)
  count <- 1
  for (i in 1:d) {
    for (j in i:d) {
      v[count] <- V[j,i]
      count <- count+1
    }
  }
  return(v)
}


vec <- function(V) {
  d <- nrow(V)
  v <- rep(0,d*d)
  count <- 1
  for (i in 1:d) {
    for (j in 1:d) {
      v[count] <- V[i,j]
      count <- count+1
    }
  }
  return(v)
}

# generate file name base
outnamebase <- function(rep) {
  o <- format(rep)
  if (rep < 10) {
    o <- paste("000",o,sep="")
  } else {
    if (rep < 100) {
      o <- paste("00",o,sep="")
    } else {
      if (rep < 1000) {
        o <- paste("0",o,sep="")
      }
    }
  }
  return(o)
}
