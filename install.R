# This is for binder, to install pkgs before we start
install.packages(c("ggplot2",
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
                   "tinytex"))

# so we can render the vignette in Binder, we need xelatex
tinytex::install_tinytex()