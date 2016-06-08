library(plspm)

# Run the customer satisfaction example form plspm

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
colnames(inner) <- rownames(inner)

# Reflective model

reflective<- matrix(
  c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
  27,6, dimnames = list(colnames(satisfaction)[1:27],colnames(inner)))

# empty formative model

formative <- matrix(0, 6, 27, dimnames = list(colnames(inner),
                                              colnames(satisfaction)[1:27]))

# Estimation using covariance matrix

satisfaction.out <- matrixpls(cov(satisfaction[,1:27]),
                           model = list(inner = inner,
                                        reflective = reflective,
                                        formative = formative))

print(satisfaction.out)


