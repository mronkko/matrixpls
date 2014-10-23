# Out of sample predictions using Model 3 from
#
# Evermann, J., & Tate, M. (2014). Out-of-Sample Predictive Ability and PLS and 
# Covariance Models. ICIS 2014 Proceedings. Retrieved from 
# http://joerg.evermann.ca/docs/EvermannTate2014b.pdf
#

if(require(simsem)){
  
  # Define the parameters to create the datasets
  
  SAMPLE <- 10000
  
  MODEL <- "
  a <~ sqrt(1/3)*a1 + sqrt(1/3)*a2 + sqrt(1/3)*a3
  b <~ sqrt(1/3)*b1 + sqrt(1/3)*b2 + sqrt(1/3)*b3
  c <~ sqrt(1/3)*c1 + sqrt(1/3)*c2 + sqrt(1/3)*c3

  k =~ 1*k1 + 1*k2 + 1*k3
  k1 ~~ .1*k1
  k2 ~~ .1*k2
  k3 ~~ .1*k3

  l =~ 1*l1 + 1*l2 + 1*l3
  l1 ~~ .1*l1
  l2 ~~ .1*l2
  l3 ~~ .1*l3

  x =~ 1*x1 + 1*x2 + 1*x3
  x1 ~~ .1*x1
  x2 ~~ .1*x2
  x3 ~~ .1*x3

  y =~ 1*y1 + 1*y2 + 1*y3
  y1 ~~ .1*y1
  y2 ~~ .1*y2
  y3 ~~ .1*y3

  z =~ 1*z1 + 1*z2 + 1*z3
  z1 ~~ .1*z1
  z2 ~~ .1*z2
  z3 ~~ .1*z3

  k ~ sqrt(.75/2)*a + sqrt(.75/2)*b
  l ~ sqrt(.75/2)*b + sqrt(.75/2)*c
  x ~ sqrt(.75)*k
  y ~ sqrt(.75/2)*k + sqrt(.75/2)*l 
  z ~ sqrt(.75)*l 

  k ~~ .25*k
  l ~~ .25*l
  x ~~ .25*x
  y ~~ .25*y
  z ~~ .25*z
  "
  
  # Generate two datasets
  
  dataSets <- sim(nRep = 2, model = MODEL, n = SAMPLE, 
                dataOnly = TRUE)
  
  
  # Estimate a PLS model with the first dataset
  
  pls.res <- matrixpls(cor(dataSets[[1]]), MODEL)
  
  # Do predictions from the second dataset
  
  predictions <- predict(pls.res, dataSets[[2]])
  
} else{
  print("This example requires the simsem package")
}

