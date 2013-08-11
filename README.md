MatrixPLS
============================

   **MatrixPLS** is a matrix based implementation of Partial Least Squares Path Modeling
   algorithm. The matrix based implementation is computationally more efficient than 
   existing PLS implementation (plspm and semPLS) and does not require raw data but
   calculates the PLS estimates from a covariance matrix. The package is designed 
   towards Monte Carlo simulations and includes functions to perform simple Monte Carlo simulations.
## Requirements and Installation

To install the stable version of **MatrixPLS** from CRAN, run in your R console:
```r
install.packages("matrixpls")
```

To install the latest development version of **MatrixPLS** from github (using the package "devtools"), run in your R console:
```
# install.packages("devtools") 
library(devtools)
install_github('matrixpls', username='mronkko')
```

TODO: Update the example

## Example Usage with a Customer Satisfaction Model 
```
library(matrixpls)

# load dataset satisfaction
data(satisfaction)

# define inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0) 
LOY = c(1,0,0,0,1,0)
sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

# define outer model list
sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

# vector of modes
sat_mod = rep("A", 6)   ## reflective indicators

# apply matrixpls with bootstrap validation
satpls = matrixpls(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)
  
# summary of results
summary(satpls)

# plot inner model results
plot(satpls, what="inner")

# plot outer model loadings
plot(satpls, what="loadings")

# plot outer model weights
plot(satpls, what="weights")
```




Author Contact
--------------
Mikko R??nkk?? (mikko.ronkko@aalto.fi)
