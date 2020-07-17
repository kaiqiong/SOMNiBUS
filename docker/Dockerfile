FROM rocker/verse:3.6.3
#Author: simon@mims.ai
#Date created: June 22th, 2020

COPY scripts/ scripts/
RUN Rscript scripts/add_tex_packages.R

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	libpoppler-cpp-dev \
	libtesseract-dev \
	libleptonica-dev \
	tesseract-ocr-eng\
	libmagick++-dev \
	libavfilter-dev \
	libwebp-dev \
	libgdal-dev \
	cargo

RUN install2.r --error \
	--deps TRUE \
	mgcv \
	testthat \
	magick \
	rhub \
	&& R -e "BiocManager::install(c('BiocStyle','BiocCheck'), update=FALSE, ask=FALSE)" \
	&& rm -rf /tmp/downloaded_packages

ENV TZ=America/Montreal
