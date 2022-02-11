# Rocker rstudio:4.1.1 Rstudio image
FROM rocker/rstudio:4.1.1
MAINTAINER Ines Assum <ines.assum@helmholtz-muenchen.de>

## Additional libraries
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libxml2-dev \
    #apt-transport-https \
    #dirmngr \
    ca-certificates \
    curl \
    gnupg2 \
    build-essential \
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
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get clean \
  && apt-get purge

RUN Rscript -e 'install.packages("devtools")' \
    && Rscript -e 'library(devtools)'

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' -e 'BiocManager::install(c("fgsea", "BiocStyle"))'\
  && Rscript -e 'install.packages(c("ranger", "ggplot2", "ggpubr", "pROC"), repos = "https://cran.rstudio.com/")'

RUN Rscript -e 'install.packages(c("formatr", "knitr", "rmarkdown", "runit", "testthat"), repos = "https://cran.rstudio.com/")'

RUN Rscript -e 'install.packages(c("caTools", "bitops", "rprojroot", "littler", "binb", "linl", "pinp", "tint"), repos = "https://cran.rstudio.com/")'

RUN Rscript -e 'install.packages(c("optparse", "XML", "plyr", "reshape2", "R.utils"), repos = "https://cran.rstudio.com/")'

RUN Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL.data/GenABEL.data_1.0.0.tar.gz", repos = NULL, type = "source")' \
  && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.8-0.tar.gz", repos = NULL, type = "source")'

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

RUN Rscript -e 'library(devtools)' -e 'devtools::install_github("andreyshabalin/MatrixEQTL", force=T)' -e 'devtools::install_github("matthiasheinig/eQTLpipeline", force=T)'

RUN Rscript -e 'installed.packages()[, "Version"]'

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    x11-apps \
    libtcl8.6 \
    libtk8.6 \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get clean \
  && apt-get purge

RUN Rscript -e 'install.packages("RandomFields", repos = "https://cran.rstudio.com/")' \
  #&&Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/mvtnorm/mvtnorm_1.0-8.tar.gz")' \
  && Rscript -e 'install.packages("Umpire", repos = "https://cran.rstudio.com/")' \
  && Rscript -e 'install.packages("Boruta", repos = "https://cran.rstudio.com/")'

RUN Rscript -e 'install.packages("geoR", repos = "https://cran.rstudio.com/")' \
  && Rscript -e 'library(devtools)' -e 'devtools::install_github("szymczak-lab/PathwayGuidedRF", force=T)' \
  && Rscript -e '.libPaths()' \
  && Rscript -e 'installed.packages()[c("ranger", "Boruta", "Umpire", "MASS", "stats"), c("Version", "LibPath")]'

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' -e 'BiocManager::install(c("msigdbr","mgsa"))'
RUN Rscript -e 'install.packages(c("readxl", "filesstrings"))'

#RUN Rscript -e 'installed.packages()[, "Version"]'

# R scripts: later imported via GitHub
RUN mkdir /data

COPY /library/* /usr/bin/
COPY /library/binaries/* /usr/bin/
COPY /data/example_data/* /data/example_data/
COPY /data/gmt/* /data/gmt/

RUN chmod a+x -R /usr/bin

EXPOSE 8787

CMD ["/init"]
