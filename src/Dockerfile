# Rocker r-ver:3.5.2 base image 
FROM rocker/r-ver:3.5.2
MAINTAINER Ines Assum <ines.assum@helmholtz-muenchen.de>


ARG RSTUDIO_VERSION
ENV RSTUDIO_VERSION=${RSTUDIO_VERSION:-1.1.463}
ARG S6_VERSION
ENV S6_VERSION=${S6_VERSION:-v1.21.7.0}
ENV S6_BEHAVIOUR_IF_STAGE2_FAILS=2
ENV PATH=/usr/lib/rstudio-server/bin:$PATH

## Download and install RStudio server & dependencies
## Attempts to get detect latest version, otherwise falls back to version given in $VER
## Symlink pandoc, pandoc-citeproc so they are available system-wide
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    file \
    git \
    libapparmor1 \
    libxml2-dev \
    libcurl4-openssl-dev \
    libedit2 \
    libssl-dev \
    lsb-release \
    psmisc \
    procps \
    python-setuptools \
    sudo \
    wget \
  && wget -O libssl1.0.0.deb http://ftp.debian.org/debian/pool/main/o/openssl/libssl1.0.0_1.0.1t-1+deb8u8_amd64.deb \
  && dpkg -i libssl1.0.0.deb \
  && rm libssl1.0.0.deb \
  && RSTUDIO_LATEST=$(wget --no-check-certificate -qO- https://s3.amazonaws.com/rstudio-server/current.ver) \
  && [ -z "$RSTUDIO_VERSION" ] && RSTUDIO_VERSION=$RSTUDIO_LATEST || true \
  && wget -q http://download2.rstudio.org/rstudio-server-${RSTUDIO_VERSION}-amd64.deb \
  && dpkg -i rstudio-server-${RSTUDIO_VERSION}-amd64.deb \
  && rm rstudio-server-*-amd64.deb \
  ## Symlink pandoc & standard pandoc templates for use system-wide
  && ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc /usr/local/bin \
  && ln -s /usr/lib/rstudio-server/bin/pandoc/pandoc-citeproc /usr/local/bin \
  && git clone https://github.com/jgm/pandoc-templates \
  && mkdir -p /opt/pandoc/templates \
  && cp -r pandoc-templates*/* /opt/pandoc/templates && rm -rf pandoc-templates* \
  && mkdir /root/.pandoc && ln -s /opt/pandoc/templates /root/.pandoc/templates \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/ \
  ## RStudio wants an /etc/R, will populate from $R_HOME/etc
  && mkdir -p /etc/R \
  ## Write config files in $R_HOME/etc
  && echo '\n\
    \n# Configure httr to perform out-of-band authentication if HTTR_LOCALHOST \
    \n# is not set since a redirect to localhost may not work depending upon \
    \n# where this Docker container is running. \
    \nif(is.na(Sys.getenv("HTTR_LOCALHOST", unset=NA))) { \
    \n  options(httr_oob_default = TRUE) \
    \n}' >> /usr/local/lib/R/etc/Rprofile.site \
  && echo "PATH=${PATH}" >> /usr/local/lib/R/etc/Renviron \
  ## Need to configure non-root user for RStudio
  && useradd rstudio \
  && echo "rstudio:rstudio" | chpasswd \
  && mkdir /home/rstudio \
  && chown rstudio:rstudio /home/rstudio \
  && addgroup rstudio staff \
  ## Prevent rstudio from deciding to use /usr/bin/R if a user apt-get installs a package
  &&  echo 'rsession-which-r=/usr/local/bin/R' >> /etc/rstudio/rserver.conf \
  ## use more robust file locking to avoid errors when using shared volumes:
  && echo 'lock-type=advisory' >> /etc/rstudio/file-locks \
  ## configure git not to request password each time
  && git config --system credential.helper 'cache --timeout=3600' \
  && git config --system push.default simple \
  ## Set up S6 init system
  && wget -P /tmp/ https://github.com/just-containers/s6-overlay/releases/download/${S6_VERSION}/s6-overlay-amd64.tar.gz \
  && tar xzf /tmp/s6-overlay-amd64.tar.gz -C / \
  && mkdir -p /etc/services.d/rstudio \
  && echo '#!/usr/bin/with-contenv bash \
          \n## load /etc/environment vars first: \
        \n for line in $( cat /etc/environment ) ; do export $line ; done \
          \n exec /usr/lib/rstudio-server/bin/rserver --server-daemonize 0' \
          > /etc/services.d/rstudio/run \
  && echo '#!/bin/bash \
          \n rstudio-server stop' \
          > /etc/services.d/rstudio/finish \
  && mkdir -p /home/rstudio/.rstudio/monitored/user-settings \
  && echo 'alwaysSaveHistory="0" \
          \nloadRData="0" \
          \nsaveAction="0"' \
          > /home/rstudio/.rstudio/monitored/user-settings/user-settings \
  && chown -R rstudio:rstudio /home/rstudio/.rstudio

RUN apt-get update \
&& apt-get install -y software-properties-common \
    apt-transport-https \
    dirmngr \
    ca-certificates \
    curl \
    gnupg2 \
    gpg \
    libc6-i386 \
    libusb-1.0-0 \
    libqt5core5a \
    libqt5serialport5 \
    build-essential \
    python-pip \
    autoconf \
    libtool \
    automake \
    cmake \
  && apt-key adv --no-tty --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF \
  && apt-get install -y mono-devel \   
  && rm -rf /var/lib/apt/lists/* \
  && apt-get clean \
  && apt-get purge

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
                ghostscript \
                lmodern \
                pandoc-citeproc \
                qpdf \
                #texlive-full \
                unzip \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get clean \
  && apt-get purge

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' -e 'BiocManager::install("fgsea", version = "3.8")' -e 'BiocManager::install("BiocStyle")'\
  && Rscript -e 'install.packages(c("ranger", "ggplot2", "ggpubr", "pROC"), repos = "https://cran.rstudio.com/")'

RUN Rscript -e 'install.packages(c("devtools", "formatr", "knitr", "rmarkdown", "runit", "testthat"), repos = "https://cran.rstudio.com/")'

RUN Rscript -e 'install.packages(c("caTools", "bitops", "rprojroot", "littler", "binb", "linl", "pinp", "tint"), repos = "https://cran.rstudio.com/")'

RUN Rscript -e 'install.packages(c("optparse","XML","plyr","reshape2","R.utils"), repos = "https://cran.rstudio.com/")'

RUN wget https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz \
  && wget https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz \
  && Rscript -e 'install.packages("GenABEL.data_1.0.0.tar.gz", repos = NULL, type = "source")' \
  && Rscript -e 'install.packages("GenABEL_1.8-0.tar.gz", repos = NULL, type = "source")' \
  && rm -rf GenABEL.data_1.0.0.tar.gz \
  && rm -rf GenABEL_1.8-0.tar.gz \
  && Rscript -e 'install.packages("MatrixEQTL", repos = "https://cran.rstudio.com/")'

RUN apt-get install -y unzip \
  && wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20181202.zip \
  && unzip plink_linux_x86_64_20181202.zip -d plink_linux_x86_64_20181202/ \
  && cd plink_linux_x86_64_20181202 \
  && mv plink /usr/bin/ \
  && mv prettify /usr/bin/ \
  && mv LICENSE /usr/bin/ \
  && cd ../ \
  && rm -rf plink_linux_x86_64_20181202 \
  && rm -rf plink_linux_x86_64_20181202.zip \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get clean \
  && apt-get purge

COPY userconf.sh /etc/cont-init.d/userconf

## running with "-e ADD=shiny" adds shiny server
COPY add_shiny.sh /etc/cont-init.d/add
COPY disable_auth_rserver.conf /etc/rstudio/disable_auth_rserver.conf
COPY pam-helper.sh /usr/lib/rstudio-server/bin/pam-helper

# R scripts: later imported via GitHub
RUN mkdir /data

COPY /library/* /usr/bin/
COPY /library/binaries/* /usr/bin/
COPY /data/example_data/* /data/example_data/
COPY /data/gmt/* /data/gmt/

RUN chmod a+x /usr/bin/MonaConsoleApp.exe \
&& chmod a+x /usr/bin/sayHello.R \
&& chmod a+x /usr/bin/enrichment_function_library.R \
&& chmod a+x /usr/bin/runner.R \
&& chmod a+x /usr/bin/multiOMICs_logregression.R \
&& chmod a+x /usr/bin/get_gmt.R

EXPOSE 8787

CMD ["/init"]
