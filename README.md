
<!-- README.md is generated from README.Rmd. Please edit that file -->

# UThwigl

[![Travis build
status](https://travis-ci.org/benmarwick/UThwigl.svg?branch=master)](https://travis-ci.org/benmarwick/UThwigl)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/benmarwick/UThwigl?branch=master&svg=true)](https://ci.appveyor.com/project/benmarwick/UThwigl)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/benmarwick/UThwigl/master?urlpath=rstudio)

The goal of UThwigl is to compute closed- and open-system
uranium-thorium (U-Th) ages of geological and archaeological samples.
Closed-system U-Th dating follows the principles presented in Edwards et
al. (2003), Richards and Dorale (2003) and Ludwig (2003). Open-system
U-Th dating is based on the diffusion-adsorption-decay (DAD) model of
Sambridge et al. (2012), which allows for advective and diffusive
transport of uranium and thorium isotopes, while including synchronous
radioactive decay.

The function, `csUTh()`, calculates closed-system ages and uses
(<sup>230</sup>Th/<sup>238</sup>U) and (<sup>234</sup>U/<sup>238</sup>U)
activity ratios on a single analysis to return an age in thousands of
years ago. The (<sup>230</sup>Th/<sup>238</sup>U) activity ratio is also
needed for detrital correction.

The function, `osUTh()`, calculates open-system ages and take the
(<sup>230</sup>Th/<sup>238</sup>U) and (<sup>234</sup>U/<sup>238</sup>U)
activity ratios collected along a transect perpendicular to the surface
of the tooth or bone and return an age in thousands of years ago.

## Web application

The package includes two web applications: one that runs `osUTh()` in
your browser where you can upload your CSV file, set the model
parameters, run the model, and view the results, and a second app that
does the same for `csUTh()`. Choose this option is you are not familiar
with R. Here’s a screenshot of the web app:

![](/Users/bmarwick/Desktop/UThwigl/docs/articles/figures/shiny-app-screenshots.png)

## Run code without downloading anything

You can run the package functions in your browser by [starting a binder
instance](https://mybinder.org/v2/gh/benmarwick/UThwigl/master?urlpath=rstudio).
This will open RStudio in your browser, together with the contents of
this GitHub repository.

## Installing the UThwigl package on your computer

To install the development version of UThwigl from GitHub on your
computer, run:

``` r
if(!require("remotes")) install.packages("remotes")
remotes::install_github("tonydoss/UThwigl")
```

Please see the [vignette](docs/articles/uthwigl.pdf) for an example of
how to use this package.

<!--
# get Tony's changes
git remote add upstream https://github.com/tonydoss/UThwigl.git
git fetch upstream
git checkout master
git merge upstream/master 

# update pkgdown
pkgdown::build_site()
-->
