# Run this once before knitting phage_datamining.Rmd
# (the Rmd also installs missing packages automatically, this script is just a convenience)

packages <- c("rmarkdown", "knitr", "readxl", "dplyr", "tidyr", "stringr", "ggplot2",
              "MASS", "pROC", "caret", "boot", "bestglm",
              "rpart", "rpart.plot", "randomForest", "gbm",
              "factoextra", "cluster", "patchwork")

missing <- packages[!packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  cat("All packages are already installed.\n")
}

