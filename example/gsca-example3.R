library(MASS)

# Run the MonteCarlo simulation from Hwang and Takane 2004

MODEL <- "
gamma1 =~ z1 + z2 + z3 + z4
gamma2 =~ z5 + z6 + z7 + z8
gamma2 ~ gamma1
"

# Data generation function

Sigma_e <- matrix(0.1,9,9)
Sigma_e[,9] <- Sigma_e[9,] <- 0
Sigma_e[5:8,1:4] <- Sigma_e[1:4,5:8] <- 0.3
diag(Sigma_e) <- 1

V <- cbind(diag(8),rep(c(0,0.3),each = 4))
W <- diag(2)[rep(c(1,2), each = 4),] * .3
A <- cbind(diag(2)[,rep(c(1,2), each = 4)] * .8,c(0.3,0))
  
Phi <- V-W %*% A

Sigma <- solve(Phi %*% t(Phi)) %*% Phi %*% Sigma_e %*% t(Phi) %*% solve(Phi %*% t(Phi))

colnames(Sigma) <- rownames(Sigma) <- paste("z",1:8,sep="")

print(cov2cor(Sigma))
stop("Sigma is weird")

generate <- function(N){
  mvrnorm(N, rep(0,8), Sigma)
}

sim1 <- matrixpls.sim(nRrep = 10, model = MODEL, n = 200,
                  outerEstimators = outer.GSCA,
                  innerEstimator = inner.GSCA,
                  generate = generate)

print(sim1)


