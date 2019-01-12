## Test environments
* local OS X install, R 3.5.1
* ubuntu 14.04 (on travis-ci), R 3.5.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* To deal with a W about ' PDF files under ‘inst/doc’', I did tools::compactPDF('vignettes/idadwigl.pdf', gs_quality = "ebook") but it makes no difference... in RStudio setting --compact-vignettes=both in build options
