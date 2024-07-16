# OBART
Ordinal Bayesian Additive Regression Trees (OBART) for ordinal outcomes.  

## Install code
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!"devtools" %in% loadedNamespaces()) {
  library(devtools)
}

if (!"remotes" %in% loadedNamespaces()) {
  library(remotes)
}

install_github("cherylwaal/OBART")

