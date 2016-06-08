# Run the education example from the book

# Sanchez, G. (2013) PLS Path Modeling with R
# Trowchez Editions. Berkeley, 2013. 
# http://www.gastonsanchez.com/PLS Path Modeling with R.pdf

education <- read.csv("http://www.gastonsanchez.com/education.csv")

Support <- c(0, 0, 0, 0, 0, 0)
Advising <- c(0, 0, 0, 0, 0, 0)
Tutoring <- c(0, 0, 0, 0, 0, 0)
Value <- c(1, 1, 1, 0, 0, 0)
# Omit two paths (compared to teh model in the book) to achieve 
# identification of the 2SLS analysis
Satisfaction <- c(0, 0, 1, 1, 0, 0)
Loyalty <- c(0, 0, 0, 0, 1, 0)

inner <- rbind(Support, Advising, Tutoring, Value, Satisfaction, Loyalty)


reflective <- diag(6)[c(rep(1,4),
                        rep(2,4),
                        rep(3,4),
                        rep(4,4),
                        rep(5,3),
                        rep(6,4)),]
formative <- matrix(0, 6, 23)

colnames(inner) <- colnames(reflective) <- rownames(formative) <- rownames(inner)
rownames(reflective) <- colnames(formative) <- colnames(education)[2:24]

model <- list(inner = inner,
              reflective = reflective,
              formative = formative)

# Reverse code two variables
education[,c("sup.under","loy.asha")] <- - education[,c("sup.under","loy.asha")]

S <- cor(education[,2:24])

# PLSc with OLS regression

plsc.res1 <- matrixpls(S,model,
                      disattenuate = TRUE,
                      parametersReflective = estimator.plscLoadings)

# PLSc with 2SLS regresssion

plsc.res2 <- matrixpls(S,model,
                      disattenuate = TRUE,
                      parametersReflective = estimator.plscLoadings,
                      parametersInner = estimator.tsls)


# Disattenuated regression with unit-weighted scales and exploratory factor analysis
# reliability estimates (with unconstrained MINRES estimator)

plsc.res3 <- matrixpls(S,model,
                       disattenuate = TRUE,
                       weightFunction = weight.fixed,
                       parametersReflective = estimator.efaLoadings)

# Disattenuated GSCA with 2SLS regression after disattenuated based on 
# confirmatory factor analysis reliability estimates


plsc.res4 <- matrixpls(S,model,
                       disattenuate = TRUE,
                       innerEstimator = inner.gsca,
                       outerEstimators = outer.gsca,
                       parametersInner = estimator.tsls,
                       parametersReflective = estimator.cfaLoadings)


# Compare the results

cbind(PLSc = plsc.res1, PLSc_2sls = plsc.res2, 
      DR = plsc.res3,.gscac = plsc.res4)

# Compare the reliability estimates

cbind(PLSc = attr(plsc.res1,"Q"), PLSc_2sls = attr(plsc.res2,"Q"), 
      DR = attr(plsc.res3,"Q"), GSCAc = attr(plsc.res4,"Q"))
