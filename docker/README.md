# unsing rocker/verse docker

https://hub.docker.com/r/rocker/verse

# to make docker
to build docker image with all necessary dependencies
'''
$ bash make_docker.sh
'''

to build docker image and push it on happyregistry
'''
$ bash make_docker.sh admin
'''

# Using the rstudio from rocker/verse container

https://hub.docker.com/r/rocker/rstudio/
bash run_docker.sh


# to run Rstudio



# to use Rstudio within your web browser

http://127.0.0.1:8787/

# once connected to RStudio

## formating .R files
from R console
configure to save file after formating
'''
> save_after_styling = TRUE
'''
format all files within a dir
'''
> tidy_dir("R/",indent = getOption("formatR.indent", 4),width.cutoff = 70)
'''
or only one file at the time
'''
> tidy_file("R/",indent = getOption("formatR.indent", 4),width.cutoff = 70)
'''

## package build
from terminal
'''
$ R CMD build .
'''

## check package tarball
from terminal
'''
$ R CMD check SOMNiBUS_0.1.0.tar.gz
'''

## BiocCheck package tarball
'''
$ R CMD BiocCheck SOMNiBUS_0.1.0.tar.gz
'''

## usefull command from devtools
from R console
load actual package
'''
> devtools::load_all(".")
'''
rendering the documentation
'''
> devtools::document(roclets = c('namespace'))
'''
run all package's tests
'''
> devtools::test()
'''

## testing
from R console
testing specific file like tests/testhat/test_groundzero.R ensuring we get same results over and over
'''
> testthat::test_file("tests/testthat/test_groundzero.R")
'''
for computing coverage
'''
> library(covr)
> report()
'''

## other usefull tricks
from R console
to shorten compilation time of vignettes that depend on internet access, it is possible to precompile the vignettes

'''
>library(knitr)
>knit("vignettes/vignette.Rmd.orig", "vignettes/vignette.Rmd")
'''
where vignette orig is the copy of the original vignette
updates are made to vignette.Rmd.orig prior to appyl precompilation
then the vignette is compiled
'''
>library(devtools)
>build_vignettes()
'''

multi platform build/check package from rhub
'''
>devtools::check_rhub()
'''

adding latex packages
'''
>tinytex::tlmgr_install(c("symbol"))
'''
