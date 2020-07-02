cores <- getOption("mc.cores")
options(mc.cores=2)

# Run the example from plspm package

# load dataset satisfaction
data(satisfaction)
# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
# outer model list
sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
# vector of modes (reflective indicators)
sat_mod = rep("A", 6)

# apply matrixpls
matrixpls.res <- matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod,
                                 scaled=FALSE, boot.val=FALSE)

print(summary(matrixpls.res))

options(mc.cores=cores)

