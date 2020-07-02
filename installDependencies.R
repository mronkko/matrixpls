# Matrixpls dependencies
install.packages(c("assertive",
"matrixcalc",
"roxygen2",
"devtools",
"psych",
"lavaan",
"simsem",
"RUnit",
"semPLS",
"boot",
"sem",
"rticles",
"stringr",
"xtable",
"R.rsp"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install("ASGSCA")
