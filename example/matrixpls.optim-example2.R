S <- diag(6)

S[upper.tri(S, diag = FALSE)] <- c(.3,
                                   -.4,-.4,
                                   .4,.4,.3,
                                   .3,.3,.3,.3,
                                   .3,.3,.3,.3,.3)

S[lower.tri(S, diag = FALSE)] <- t(S)[lower.tri(S, diag = FALSE)]

inner <- matrix(c(0,0,1,
                  0,0,1,
                  0,0,0),3,3)

reflective <- diag(3)[rep(1:3, each = 2),]

formative <- matrix(0,3,6)

colnames(inner) <- rownames(inner) <- colnames(reflective) <-
  rownames(formative) <- c("A","B","C")

colnames(S) <- rownames(S) <- colnames(formative) <-
  rownames(reflective) <- c("a1","a2","b1","b2","c1","c2")


model <- list(inner = inner,
              reflective = reflective,
              formative = formative)

modeA <- matrixpls(S,model, outerEstim = outerEstim.modeA)

modeB <- matrixpls(S,model, outerEstim = outerEstim.modeB)

fixed <- matrixpls(S,model, weightFun = weightFun.fixed)

optimR2 <- matrixpls(S,model, weightFun = weightFun.optim)

optimGoF <- matrixpls(S,model, weightFun = weightFun.optim,
                      optimCrit = function(matrixpls.res){
                        - gof(matrixpls.res)
                      })


# R2
rbind(ModeA = r2(modeA),
      ModeB = r2(modeB),
      Fixed = r2(fixed),
      Optim = r2(optimR2))

# GoF abs

rbind(ModeA = gof(modeA),
      ModeB = gof(modeB),
      Fixed = gof(fixed),
      Optim = gof(optimGoF))
