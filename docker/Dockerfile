FROM rocker/verse:latest
#Author: simon@mims.ai
#Date created: June 22th, 2020

RUN install2.r --error \
    --deps TRUE \
    mgcv \
    testthat\
    && R -e "BiocManager::install(c('BiocStyle','BiocCheck'), update=FALSE, ask=FALSE)"

ENV TZ=America/Montreal
