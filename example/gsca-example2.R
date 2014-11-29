# Run the example from http://www.sem-gesca.org/example.php

if(!exists("maleData")){
  data <- read.table("http://www.sem-gesca.org/rick2_0.dat", header = TRUE)
  maleData <- data[data$gender == 1,-1]
}

inner <- matrix(0,4,4)
inner[c(2,7,8)] <- 1

reflective <- matrix(0,21,4)
reflective[1:8,1] <- 1
reflective[9:14,2] <- 1
reflective[15:18,3] <- 1
reflective[19:21,4] <- 1

formative <- matrix(0,4,21)

colnames(formative) <- rownames(reflective) <- names(maleData)
colnames(inner) <- rownames(inner) <- rownames(formative) <- colnames(reflective) <- 
  c("Organizational prestige", "Organizational indentification",
    "Affective commitment (joy)", "Affective commitment (love)")

model <- list(inner = inner, 
              reflective = reflective,
              formative = formative)

# Estimate using alternating least squares

matrixpls.res2 <- matrixpls(cov(maleData),  model,
                            outerEstimators = outer.GSCA,
                            innerEstimator = inner.GSCA)


# Estimate using direct minimization of the estimation criterion

matrixpls.res1 <- matrixpls(cov(maleData),  model,
                            weightFunction = weight.optim,
                            optimCriterion = optim.GSCA)


# Compare the weights
print(t(attr(matrixpls.res1,"W")))
print(t(attr(matrixpls.res2,"W")))

print(t(attr(matrixpls.res1,"W"))-t(attr(matrixpls.res2,"W")))

# Compare the paths
print(t(attr(matrixpls.res1,"beta")))
print(t(attr(matrixpls.res2,"beta")))

print(t(attr(matrixpls.res1,"beta"))-t(attr(matrixpls.res2,"beta")))


# Compare the loadings
print(loadings(matrixpls.res1))
print(loadings(matrixpls.res2))


print(loadings(matrixpls.res1)-loadings(matrixpls.res2))




