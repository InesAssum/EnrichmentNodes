# Ubuntu base image 
FROM ubuntu:18.04
MAINTAINER Ines Assum <ines.assum@helmholtz-muenchen.de>

# Setting variables to prevent errors at R installation
ENV TZ=Europe/Minsk
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

#installation of software
RUN apt-get update && apt-get install -y software-properties-common \
&& add-apt-repository ppa:git-core/ppa \
&& add-apt-repository ppa:ubuntu-toolchain-r/test \
&& apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF \
&& echo "deb https://download.mono-project.com/repo/ubuntu stable-bionic main" | tee /etc/apt/sources.list.d/mono-official-stable.list \
&& apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
&& apt update \
&& apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    mercurial \
    curl \
    wget \
    libcurl4-openssl-dev \
    pkg-config \
    libfreetype6-dev \
    libmysqlclient-dev \
    zlib1g-dev \
    tcsh \
    cmake \
    r-base \
    libbz2-dev \
    liblzma-dev \
    liblzo2-dev \
    locales \
    ca-certificates \
    mono-complete \
    monodevelop \
#&& update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc 60 --slave /usr/bin/g++ g++ /usr/bin/g++ \
&& rm -rf /var/lib/apt/lists/* \
&& apt-get clean \
&& apt-get purge


# Here a few code chunks on how to install R packages
# CRAN packages:
# RUN Rscript -e "install.packages('PACKAGENAME')"

# Bioconductor packages - until R version 3.4?
# RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("fgsea")'
#
# Bioconductor packages - starting R version 3.5
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("fgsea", version = "3.8")
# RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' -e 'BiocManager::install("fgsea", version = "3.8")'


RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("fgsea")'
RUN Rscript -e 'install.packages("ggplot2")'
RUN Rscript -e 'install.packages("ggpubr")'

RUN echo '#!/usr/bin/env Rscript \n\
args = commandArgs(trailingOnly=TRUE) \n\
if (!(length(args)==2)) { \n\
  stop("One argument must be supplied ((input file).n)string)", call.=FALSE) \n\
} \n\
sayHello <- function(){ \n\
  writeLines(args[1], args[2]) \n\
} \n\
sayHello()' > /usr/bin/sayHello.R# Ubuntu base image 
FROM ubuntu:18.04
MAINTAINER Ines Assum <ines.assum@helmholtz-muenchen.de>

# Setting variables to prevent errors at R installation
ENV TZ=Europe/Minsk
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

#installation of software
RUN apt-get update && apt-get install -y software-properties-common \
&& apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF \
&& echo "deb https://download.mono-project.com/repo/ubuntu stable-bionic main" | tee /etc/apt/sources.list.d/mono-official-stable.list \
&& apt update \
&& apt-get update \
&& apt-get install -y \
    r-base \
    monodevelop \
#   mono-complete \    
&& rm -rf /var/lib/apt/lists/* \
&& apt-get clean \
&& apt-get purge

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("fgsea")'
#RUN Rscript -e 'install.packages("ggplot2")'
#RUN Rscript -e 'install.packages("ggpubr")'

# R scripts: later imported via GitHub

RUN mkdir /data

COPY /files/* /usr/bin/
COPY /files/Release/* /usr/bin/
COPY /data/example_data/* /data/example_data/
COPY /data/gmt/* /data/gmt/

RUN chmod a+x /usr/bin/MonaConsoleApp.exe \
&& chmod a+x /usr/bin/sayHello.R \
&& chmod a+x /usr/bin/run_fgsea.R \
&& chmod a+x /usr/bin/run_mona_single.R
