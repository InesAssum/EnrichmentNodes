# Ubuntu base image 
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
&& rm -rf /var/lib/apt/lists/* \
&& apt-get clean \
&& apt-get purge

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("fgsea")'
#RUN Rscript -e 'install.packages("ggplot2")'
#RUN Rscript -e 'install.packages("ggpubr")'

# R script: sayHello.R, later imported via GitHub
RUN echo '#!/usr/bin/env Rscript \n\
args = commandArgs(trailingOnly=TRUE) \n\
if (!(length(args)==2)) { \n\
  stop("One argument must be supplied ((input file).n)string)", call.=FALSE) \n\
} \n\
sayHello <- function(){ \n\
  writeLines(args[1], args[2]) \n\
} \n\
sayHello()' > /usr/bin/sayHello.R
RUN chmod a+x /usr/bin/sayHello.R

# R script: run_fgsea.R, later imported via GitHub
RUN echo '#!/usr/bin/env Rscript \n\
args = commandArgs(trailingOnly=TRUE) \n\
if (!(length(args)==2)) { \n\
  stop("One argument must be supplied ((input file).n)string)", call.=FALSE) \n\
} \n\
sayHello <- function(){ \n\
  writeLines(args[1], args[2]) \n\
} \n\
sayHello()' > /usr/bin/sayHello.R \
&& chmod a+x /usr/bin/sayHello.R

RUN echo '#!/usr/bin/env Rscript \n\
args = commandArgs(trailingOnly=TRUE) \n\
 \n\
data <- args[1] \n\
result <- args[2] \n\
gmt <- args[3] \n\
minSize <- args[4] \n\
maxSize <- args[5] \n\
nperm <- args[6] \n\
mygmt <- args[7] \n\
 \n\
if(gmt=="KEGG"){gmt <- "/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"} \n\
if(gmt=="GO"){gmt <- "/data/gmt/c5.bp.v6.0.symbols.gmt"} \n\
if(gmt=="custom"){gmt <- mygmt} \n\
 \n\
# load libraries \n\
library(fgsea) \n\
 \n\
# load data \n\
if(grep("\\.RDS | \\.rds", data)){ \n\
  rank <- readRDS(data) \n\
} else if (grep("\\.RData | \\.Rdata | \\.rdata", data)){ \n\
  # # need to test this further, incomplete (not working) \n\
  # input.env <- environment() \n\
  # load(data, envir = input.env) \n\
  # input.env$data <- ls(envir = input.env) # (not working) \n\
  # rank <- input.env$data \n\
} \n\
 \n\
pathways <- gmtPathways(gmt) \n\
GSEA <- fgsea(pathways = pathways, \n\
              stats = rank, \n\
              minSize=minSize, \n\
              maxSize=maxSize, \n\
              nperm=nperm) \n\
table[, c("Pathway", "score", "score2", "pvalue", "padjust")] <- \n\
  c(GSEA[, c("Pathway", "ES")], NA, GSEA[, c("pvalue", "padjust")]) \n\
 \n\
res <- list("summary"=table, \n\
            "res"=GSEA, \n\
            "method"="GSEA", \n\
            "opts"=list("data"=data, \n\
                        "results"=result, \n\
                        "gmt"=gmt, \n\
                        "minSize"=minSize, \n\
                        "maxSize"=maxSize, \n\
                        "nperm"=nperm)) \n\
 \n\
saveRDS(res, file=result) \n\
q(save = "no")' > /usr/bin/run_fgsea.R
RUN chmod a+x /usr/bin/run_fgsea.R


# R script: run_mona1.R, later imported via GitHub
RUN echo '#!/usr/bin/env Rscript \n\
args = commandArgs(trailingOnly=TRUE) \n\
 \n\
data <- args[1] \n\
result <- args[2] \n\
gmt2 <- args[3] \n\
minSize <- args[4] \n\
maxSize <- args[5] \n\
cutoff <- args[6] \n\
sign <- args[7] \n\
rev <- args[8] \n\
mygmt <- args[9] \n\
 \n\
if(gmt=="KEGG"){gmt <- "/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"} \n\
if(gmt=="GO"){gmt <- "/data/gmt/c5.bp.v6.0.symbols.gmt"} \n\
if(gmt=="custom"){gmt <- mygmt} \n\
 \n\
p.mono <- "mono" \n\
p.mona <- "/usr/bin/MonaConsoleApp.exe" \n\
 \n\
library(fgsea) \n\
 \n\
gmt <- gmtPathways(gmt2) \n\
hidden <- unique(unlist(gmt)) \n\
Ass <- matrix(NA, dimnames = list(hidden, names(gmt)), nrow = length(hidden), ncol = length(gmt)) \n\
for (i in 1:dim(Ass)[2]){ \n\
  Ass[,i] <- as.numeric(hidden %in% gmt[[i]]) \n\
} \n\
 \n\
if(grep("\\.RDS | \\.rds", data)){ \n\
  values <- readRDS(data) \n\
} else if (grep("\\.RData | \\.Rdata | \\.rdata", data)){ \n\
  # # need to test this further, incomplete (not working) \n\
  # input.env <- environment() \n\
  # load(data, envir = input.env) \n\
  # input.env$data <- ls(envir = input.env) # (not working) \n\
  # rank <- input.env$data \n\
} \n\
 \n\
genes <- unique(names(rank)) \n\
hidden2 <- intersect(rownames(Ass[rowSums(Ass[genes, ]) > minSize,]), \n\
                     rownames(Ass[rowSums(Ass[genes, ]) < maxSize,])) \n\
Ass2 <- Ass[hidden2, colnames(Ass)[which(colSums(Ass[hidden2,])>minSize & \n\
                                         colSums(Ass[hidden2,])<maxSize)]] \n\
 \n\
mona1 <- tempdir() \n\
p.out <- paste0(mona1, "/", "output.txt") \n\
p.assign <- paste0(mona1, "/", "assignmentMatrix.txt") \n\
p.terms <- paste0(mona1, "/", "terms.txt") \n\
p.yn <- paste0(mona1, "/", "values.txt") \n\
 \n\
if(!sign){ \n\
  assign <- apply(Ass2,1,function(x) paste(which(x>0)-1,collapse=",")) \n\
  if(rev){ \n\
    yn <- as.numeric(abs(values) > cutoff) \n\
  }else{ \n\
    yn <- as.numeric(abs(values) < cutoff) \n\
  } \n\
} else if (sign){ \n\
  A <- Ass2 \n\
  rownames(A) <- paste0(rownames(Ass2), "_up") \n\
  colnames(A) <- paste0(colnames(Ass2), "_up") \n\
  NA1 <- matrix(0, dim(Ass2)[1], dim(Ass2)[2]) \n\
  rownames(NA1) <- paste0(rownames(Ass2), "_up") \n\
  colnames(NA1) <- paste0(colnames(Ass2), "_down") \n\
  NA2 <- matrix(0, dim(Ass2)[1], dim(Ass2)[2]) \n\
  rownames(NA2) <- paste0(rownames(Ass2), "_down") \n\
  colnames(NA2) <- paste0(colnames(Ass2), "_up") \n\
  B <- Ass2 \n\
  rownames(B) <- paste0(rownames(Ass2), "_down") \n\
  colnames(B) <- paste0(colnames(Ass2), "_down") \n\
  Ass_pm <- rbind(cbind(A, NA1), cbind(NA2, B)) \n\
  assign <- apply(Ass_pm,1,function(x) paste(which(x>0)-1,collapse=",")) \n\
  terms <- colnames(Ass_pm) \n\
  if(rev){ \n\
    yn_p <- as.numeric(abs(values) > cutoff & values > 0) \n\
    yn_m <- as.numeric(abs(values) > cutoff & values < 0) \n\
    yn <- c(yn_p, yn_m) \n\
  }else{ \n\
    yn_p <- as.numeric(abs(values) < cutoff & values > 0) \n\
    yn_m <- as.numeric(abs(values) < cutoff & values < 0) \n\
    yn <- c(yn_p, yn_m) \n\
  } \n\
} \n\
 \n\
write(assign, file=p.assign) \n\
write.table(terms, file=p.terms, col.names=F, row.names = F, quote=F) \n\
write.table(yn, file=p.yn, col.names=F, row.names = F, quote=F) \n\
 \n\
tries <- 0 \n\
while (!file.exists(p.out) | tries<10){ \n\
  sys1 <- system(paste(p.mono, p.mona, "0", \n\
                       p.assign, p.yn, p.terms, p.out, "1", \n\
                       sep = " "), intern = T) \n\
  tries <- tries+1 \n\
} \n\
 \n\
res1 <- read.table(file=p.out, sep="\t", h=F) \n\
res1 <- res1[order(res1$V2, decreasing=T),] \n\
colnames(res1) <- c("pathway", "posterior") \n\
unlink(mona1, recursive = T) \n\
table[, c("Pathway", "score", "score2", "pvalue", "padjust")] <- \n\
  c(res1[, c("pathway", "posterior")], NA, NA, NA) \n\
 \n\
res <- list("summary"=table, \n\
            "res"=res1, \n\
            "method"="MONA", \n\
            "opts"=list("data"=data, \n\
                        "results"=result, \n\
                        "gmt"=gmt, \n\
                        "minSize"=minSize, \n\
                        "maxSize"=maxSize, \n\
                        "cutoff"=cutoff, \n\
                        "sign"=sign, \n\
                        "rev"=rev)) \n\
 \n\
saveRDS(res, file=result) \n\
q(save = "no")' > /usr/bin/run_mona1.R
RUN chmod a+x /usr/bin/run_mona1.R

RUN mkdir /data \
&& mkdir /data/example_data \
&& mkdir /data/gmt \
&& chmod -R a+x data

#COPY MonaConsoleApp.exe /usr/bin/MonaConsoleApp.exe \
#&& chmod -R a+x /usr/bin/MonaConsoleApp.exe




