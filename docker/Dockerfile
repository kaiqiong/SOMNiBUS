FROM rocker/verse:latest
#Author: simon@mims.ai
#Date created: June 22th, 2020

COPY scripts/ scripts/

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	libpoppler-cpp-dev \
	libtesseract-dev \
	libleptonica-dev \
	tesseract-ocr-eng\
	libmagick++-dev \
	libavfilter-dev
RUN install2.r --error \
	--deps TRUE \
	mgcv \
	testthat \
	lintr \
	styler \
	magick \
	rhub \
	&& R -e "BiocManager::install(c('BiocStyle','BiocCheck'), update=FALSE, ask=FALSE)" \
	&& rm -rf /tmp/downloaded_packages


RUN Rscript scripts/add_tex_packages.R

ENV TZ=America/Montreal
