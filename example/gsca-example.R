# Run the example from plspm package using GSCA estimation

library(plspm)

# load dataset satisfaction
data(satisfaction)

# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)

inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

# reflective relationships
reflective <- matrix(0, 27, 6)
for(lv in 1:6){
	reflective[list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)[[lv]],lv] <-1
}
# no formative relationships
formative <- matrix(0, 6, 27)

model <- list(inner=inner, reflective=reflective, formative=formative)

# apply matrixpls with GSCA estimators

matrixpls(cor(satisfaction[,1:27]), model, 
					outerEstimators = outer.GSCA, 
					innerEstimator = inner.GSCA)
