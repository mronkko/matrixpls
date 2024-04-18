matrixpls
============================

   **matrixpls** is a matrix based implementation of Partial Least Squares Path Modeling
   algorithm. The matrix based implementation is computationally more efficient than 
   existing PLS implementation (plspm and semPLS) and does not require raw data but
   calculates the PLS estimates from a covariance matrix. The package is designed 
   towards Monte Carlo simulations with simsem.
   
## Development

matrixpls has three development targets or main requirements

1. matrixpls should be able to estimate all relevant [vignettes from simsem](https://github.com/simsem/simsem/wiki/Vignette)
2. matrixpls works as drop-in replacement for the most recent version of [plspm](http://cran.r-project.org/web/packages/plspm/)
3. matrixpls works as drop-in replacement for the most recent version of [plsSEM](http://cran.r-project.org/web/packages/semPLS/)

All simsem vignettes should work by just replacing the call to sim with a call to matrixpls.sim. Compatibility with plspm and semPLS should be as easy as substituting the calls to plspm function with matrixpls.plspm and semPLS with matrixpls.semPLS.

matrixpls uses test driven development with [Runit](http://cran.r-project.org/web/packages/RUnit/). This means that the tests are written first and then software is implemented to conform with the tests.

## Requirements and Installation

To install the latest development version of **matrixpls** from github (using the package "devtools"), run in your R console:
```
# install.packages("devtools") 
library(devtools)
install_github("mronkko/matrixpls")
```

To check out the source code and install from local copy, run the following shell commands

```
git clone https://github.com/mronkko/matrixpls.git
R CMD BUILD matrixpls
R CMD CHECK matrixpls_*.tar.gz
R CMD INSTALL matrixpls_*.tar.gz
```

You may need to install development versions of Lavaan and simsem to run the development version of
matrixpls

```
 install.packages("lavaan", repos="http://www.da.ugent.be", type="source")

 install.packages("simsem", repos="http://rweb.quant.ku.edu/kran", type="source")
```
 
Author Contact
--------------
Mikko Rönkkö (mikko.ronkko@jyu.fi)
