# These commands can be used to chekc the package before submitting to CRAN

devtools::check(vignettes = FALSE, env_vars=c(`_R_CHECK_DEPENDS_ONLY_` = "true")) 