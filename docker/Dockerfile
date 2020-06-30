FROM rocker/tidyverse:latest
#Author: simon@mims.ai
#Date created: June 22th, 2020

RUN install2.r mgcv pacman testthat BiocManager && R -e "BiocManager::install(c('BiocStyle','BiocCheck'), update=FALSE, ask=FALSE)"

ENV TZ=America/Montreal
