# Matrixpls dependencies
install.packages(c("assertive",
"roxygen2",
"devtools",
"psych",
"lavaan",
"simsem",
"RUnit",
# "semPLS",
"boot",
"sem",
"rticles",
"stringr",
"xtable",
"R.rsp",
"qpdf"))

# Some documentation is inherited from the abandoned semPLS package

devtools::install_url("https://cran.r-project.org/src/contrib/Archive/semPLS/semPLS_1.0-10.tar.gz")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BioConductor does not always work with R-devel. Try using BioConductor and
# fall back to git install if this fails.

tryCatch(
  {
    BiocManager::install("ASGSCA")
  },
  error=function(cond) {
    message(cond)
    install_bioc("ASGSCA")
  }
) 
