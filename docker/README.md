# unsing rocker/verse docker

https://hub.docker.com/r/rocker/verse

# to make docker

bash make_docker.sh

# to run Rstudio

bash run_docker.sh

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
