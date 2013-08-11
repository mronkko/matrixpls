library('RUnit')
 
source('R/matrixpls.sim.R')
 
test.suite <- defineTestSuite("matrixpls",
                              dirs = file.path("tests"),
                              testFileRegexp = '^\\d+\\.R')
 
test.result <- runTestSuite(test.suite)
 
printTextProtocol(test.result)