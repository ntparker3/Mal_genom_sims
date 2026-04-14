FROM ubuntu:24.04
LABEL authors="Alfred Hubbard"

ARG DEBIAN_FRONTEND="noninteractive"
ARG SHELL="/bin/bash"
ARG LANG="en_US.UTF-8"
ARG LANGUAGE="en_US.UTF-8"
ARG LC_ALL="en_US.UTF-8"
ARG CPU_COUNT=5
ARG TIME_ZONE=Etc/UTC

RUN ulimit -n 10000

# install basics
RUN apt-get update
RUN apt-get -yq dist-upgrade
RUN apt-get install -yq --no-install-recommends --fix-missing \
    autotools-dev \
    autoconf \
    libtool \
    automake \
    build-essential \
    curl \
    wget \
    file \
    git \
    locales \
    libssl-dev \
    libcurl4-gnutls-dev \
    ca-certificates \
    xz-utils \
    zlib1g-dev \
    libbz2-dev \
    liblzma5 \
    liblzma-doc \
    liblzma-dev \
    openssh-client \
    python3-dev \
    python3-pip \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libmpfr-dev \
    libgmp3-dev \
    libudunits2-dev \
    cmake \
    gdal-bin \
    libgdal-dev \
    pandoc \
    fontconfig \
    fonts-dejavu \
    fonts-freefont-ttf

# Configure locale
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen && \
    update-locale LANG=en_US.UTF-8

# Set environment variables
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8  
ENV LC_ALL=en_US.UTF-8

# Rebuild font cache
RUN fc-cache -fv

# install R
RUN apt-get install -yq --no-install-recommends r-base r-base-dev

RUN DEBIAN_FRONTEND=noninteractive TZ=$TIME_ZONE apt-get -y install tzdata

# R configuration
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = ${CPU_COUNT})" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages(c("remotes"))'
## attempt to load libraries to make sure they installed
RUN R -e 'library("remotes")'

# Basic R packages
RUN Rscript -e "remotes::install_cran(c(\
    'dplyr', \
    'forcats', \
    'ggplot2', \
    'knitr', \
    'lubridate', \
    'magrittr', \
    'purrr', \
    'readr', \
    'readxl', \
    'stringr', \
    'tibble', \
    'tidyr', \
    'rmarkdown', \
    'optparse', \
    'patchwork', \
    'doParallel', \
    'parallelly'\
    ), Ncpus = ${CPU_COUNT})"

# attempt to load libraries to make sure they installed
RUN R -e 'library("dplyr")'
RUN R -e 'library("forcats")'
RUN R -e 'library("ggplot2")'
RUN R -e 'library("knitr")'
RUN R -e 'library("lubridate")'
RUN R -e 'library("magrittr")'
RUN R -e 'library("purrr")'
RUN R -e 'library("readr")'
RUN R -e 'library("readxl")'
RUN R -e 'library("stringr")'
RUN R -e 'library("tibble")'
RUN R -e 'library("tidyr")'
RUN R -e 'library("rmarkdown")'
RUN R -e 'library("optparse")'
RUN R -e 'library("patchwork")'
RUN R -e 'library("doParallel")'
RUN R -e 'library("parallelly")'

# Additional R packages
# For some reason I get an error installing the most recent version of 
# Dcifer - I still need to dig in and troubleshoot
RUN R -e "remotes::install_github('EPPIcenter/dcifer@v1.3.1')"
RUN R -e 'library("dcifer")'
RUN Rscript -e "remotes::install_cran('argparse')"
RUN R -e 'library("argparse")'
RUN Rscript -e "remotes::install_github('PlasmoGenEpi/recombuddy@develop')"
RUN R -e 'library("recombuddy")'
