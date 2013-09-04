library(RUnit)
library(parallel)
library(matrixpls)

options(boot.ncpus = detectCores())
options(boot.parallel = "multicore")

 
test.suite <- defineTestSuite("matrixpls",
                              dirs = file.path("tests"),
                              testFileRegexp = '^\\d+\\.R')
 
test.result <- runTestSuite(test.suite)
 
printTextProtocol(test.result)