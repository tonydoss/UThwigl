# This is for binder, to install pkgs before we start
install.packages(c("ggplot2",
                   "deSolve",
                   "cowplot",
                   "knitr",
                   "rmarkdown",
                   "devtools",
                   "dplyr",
                   "bookdown",
                   "xtable",
                   "here",
                   "rticles",
                   "testthat",
                   "tinytex",
                   "roxygen2"))

# so we can render the vignette in Binder, we need xelatex
tinytex::install_tinytex()

