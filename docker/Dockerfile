FROM rocker/verse:4.0.2
#Author: simon@mims.ai
#Date created: June 22th, 2020

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
	cargo \
	&& tlmgr install \
	ec \
	graphics \
	symbol \
	&& install2.r --error \
	--deps TRUE \
	mgcv \
	VGAM \
	testthat \
	magick \
	rhub \
	&& R -e "BiocManager::install(c('BiocStyle','BiocCheck'), update=FALSE, ask=FALSE)" \
	&& rm -rf /tmp/downloaded_packages

ENV TZ=America/Montreal
