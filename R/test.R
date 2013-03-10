#
# Creates Monte Carlo study with four latent variables, A-D each with 3 indicators and estimates
#
# Set up parallel bootstrapping

options(boot.parallel="multicore")
options(boot.ncpus=detectCores())

# Create a Monte Carlo study

Phi<-diag(c(1,1, #Variances of A and B
			1-(.5*(.5+.3*.1)+.1*(.1+.3*.5)), #Error variance of C
			1-2*.3*(.5+.3*.1), # Error variance of D
			rep(1-.7^2,12))) # Indicator error variances
			
Phi[1,2]<-Phi[2,1]<-.3 #Correlation between A and B

colnames(Phi)<-c("A","B","eC","eD",paste("ex",1:12,sep=""))
rownames(Phi)<-colnames(Phi)

gamma <- matrix(0,14,16)
gamma[1,1] <- .5 # C on A
gamma[1,2] <- .1 # C on B
gamma[2,1] <- .3 # D on A
gamma[3:5,1] <- .7 # x1-x3 on A
gamma[6:8,2] <- .7 # x4-x6 on B

gamma[1:14,3:16] <- diag(14) # errors

colnames(gamma)<-colnames(Phi)
rownames(gamma)<-c("C","D",paste("x",1:12,sep=""))

beta <- matrix(0,14,14)
beta[2,1] <- .3 # D on C
beta[9:11,1] <- .7 # x7-x9 on C
beta[12:14,2] <- .7 # x10-x12 on D

rownames(beta)<-rownames(gamma)
colnames(beta)<-rownames(beta)

inner<-matrix(0,4,4)
inner[c(3,4),1]<-1
inner[c(3,4),2]<-1
inner[4,3]<-1

matrixpls.mc(rep=10, beta=beta, gamma=gamma, Phi=Phi, observed.dependents=c(3:14), n=100, plspm = TRUE, 
	inner = inner, outer = split(1:12,rep(1:4, each=3)), modes = rep("A",4), boot.val=TRUE, br=500)


