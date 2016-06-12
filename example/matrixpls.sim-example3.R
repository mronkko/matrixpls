library(simsem)

# Specify the data generation model using SimSem syntax

# Loadings
LY <- bind(matrix(c(.6, .7, .8, 0, 0, 0,
                  0, 0, 0, .6, .7, .8), 6,2))

# Factor correlations (2 different values)
RPS1 <- binds(matrix(c(1, 0, 0, 1), 2, 2))
RPS2 <- binds(matrix(c(1, .3, .3, 1), 2, 2))

# Error terms
RTE <- binds(diag(c(.64, .51, .36, .64, .51, .36)))

generateModels <- list(model.cfa(LY = LY, RPS = RPS1, RTE = RTE),
                       model.cfa(LY = LY, RPS = RPS2, RTE = RTE))

# Specify the computed model with lavaan syntax
analyzeModel <- "f1=~ y1 + y2 + y3
                 f2=~ y4 + y5 + y6
                 f2~f1"


coefficients <- list()
reliabilities <- list() 

# We have 3 weight techniques and 2 populations: 6 conditions

for(i in 1:6){
  population <- rep(1:2, each=3)[i]
  technique <-  rep(1:3, 2)[i]
  out <- matrixpls.sim(nRep = 500, model = analyzeModel,
                       weightFunction = ifelse(technique == 1, weight.fixed, weight.pls),
                       # Use Mode B for the third estimator, otherwise use Mode A
                       # This argument is ignored by fixed weights
                       outerEstimators = ifelse(technique == 3, outer.modeB, outer.modeA),
                       n = 100, generate = generateModels[[population]],
                       sequential = TRUE, saveLatentVar = TRUE,
                       multicore = TRUE, boot.R = FALSE, fitIndices = NULL)

  coefficients[[i]] <- out@coef[,1]
  # Extract the reliability of the first composite from the matrixpls 
  # results objects stored as extra output
  reliabilities[[i]] <- unlist(lapply(out@extraOut, function(x) attr(x,"R")[1]))
}

par(mfrow = c(2,2))
par(mar = c(4,4,1,0.5))

plot(density(coefficients[[1]]), main = "Beta = 0", xlab = "Estimate", 
     xlim = c(-.6,.6), ylim = c(0,6))
lines(density(coefficients[[2]]), col = "blue")
lines(density(coefficients[[3]]), col = "red")
legend("topleft", c("Unit weights", "PLS Mode A", "PLS Mode B"), 
       col=c("black","blue","red"), lty = 1, cex=0.5)

plot(density(coefficients[[4]]), main = "Beta = .3", xlab = "Estimate", 
     xlim = c(-.6,.6), ylim = c(0,6))
lines(density(coefficients[[5]]), col = "blue")
lines(density(coefficients[[6]]), col = "red")

plot(density(reliabilities[[1]]), main = "", xlab = "Reliability",
     xlim = c(-1,1), ylim = c(0,9))
lines(density(reliabilities[[2]]), col = "blue")
lines(density(reliabilities[[3]]), col = "red")

plot(density(reliabilities[[4]]), main = "", xlab = "Reliability",
     xlim = c(-1,1), ylim = c(0,9))
lines(density(reliabilities[[5]]), col = "blue")
lines(density(reliabilities[[6]]), col = "red")
