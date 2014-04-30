# Runs the second part of
#
# Figure 2. Distribution of reliability for partial least squares (PLS) Mode A, PLS Mode B, and 
# summed scales in the two-construct model over 500 replications.
#
# from
# 
# Rönkkö, M., & Evermann, J. (2013). A Critical Examination of Common Beliefs About Partial Least 
# Squares Path Modeling. Organizational Research Methods, 16(3), 425–448. 
# doi:10.1177/1094428112474693
# 
# assessing the reliability of construct scores from PLS Mode A, PLS Mode B, and summed scales

if(require(simsem)){
  
  #Define Simulation Parameters
  
  SAMPLE <- 100
  REPLICATIONS <- 10
  
  MODEL <- "
A =~ .8*a1 + .7*a2 + .6*a3
a1~~(1-.8^2)*a1
a2~~(1-.7^2)*a2
a3~~(1-.6^2)*a3
A~~1*A
B =~ .8*b1 + .7*b2 + .6*b3
b1~~(1-.8^2)*b1
b2~~(1-.7^2)*b2
b3~~(1-.6^2)*b3
B ~ .3*A
B~~(1-.3^2)*B
"
  
  modeA <- matrixpls.sim(nRep = REPLICATIONS, model = MODEL, n = SAMPLE, 
                         # Generate data by first generating exogenous variables from normal 
                         # distribution
                         sequential = TRUE,
                         # Control parameters
                         multicore = FALSE, boot.R = FALSE, stopOnError = TRUE, fitIndices = NULL)
  
} else{
  print("This example requires the simsem package")
}

