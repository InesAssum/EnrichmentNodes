#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ------------------------------------------------------------------------------
#' Function library for enrichment analysis
#'
#' @author Kristof Gilicze
#' @param data string // file path to .RDS input file: named numerical vector, names=gene symbols#
#' @param result string // file path to save results as .RDS 
#' @param gmt string // Dir containing Gene set definition to be used (KEGG/GO/custom)
#' @param minSize integer // min size of GSs to be considered, default: 15
#' @param maxSize integer // max size of GSs to be considered, default: 500
#' @param nperm integer // number of permutations, default: 10000
#' 

# ------------------------------------------------------------------------------
library("fgsea")

# ------------------------------------------------------------------------------
########################## Functions ##############################
#' Mona 1D
run_mona1 <- function(input, gmt, minSize, maxSize, cutoff, sign, rev, debug=FALSE) {
  library(fgsea)
  p.mona <- "/usr/bin/MonaConsoleApp.exe"
  p.mono <- "mono"
  # gmt prep
  gmt <- gmtPathways(gmt)
  sign <- sign=="yes"
  rev <- rev=="small"
  # gmt <- gmtPathways(mygmt)
  
  hidden <- unique(unlist(gmt))
  Ass <- matrix(NA, dimnames = list(hidden, names(gmt)), nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(Ass)[2]){
    Ass[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  # data reading
  if(grepl("\\.RDS", data) | grepl("\\.rds", data)){
    values <- readRDS(data)
  } else if (grepl("\\.RData", data) | grepl("\\.Rdata", data) | grepl("\\.rdata", data) | grepl("\\.rda", data)){}
  # filtering
  hidden2 <- unique(intersect(names(values[!is.na(values)]), hidden))
  #length(hidden2)
  Ass2 <- Ass[hidden2, colnames(Ass)[which(colSums(Ass[hidden2,]) > minSize &
                                             colSums(Ass[hidden2,]) < maxSize)]]
  # temp creation
  mona1 <- tempdir()
  dir.create(mona1)
  if(debug) print(mona1)
  p.out <- paste0(mona1, "/", "output.txt")
  p.assign <- paste0(mona1, "/", "assignmentMatrix.txt")
  p.terms <- paste0(mona1, "/", "terms.txt")
  p.yn <- paste0(mona1, "/", "values.txt")
  
  # transform and filter
  if(!sign){
    assign <- apply(Ass2,1,function(x) paste(which(x>0)-1,collapse=","))
    if(rev){
      yn <- as.numeric(abs(values[hidden2]) > cutoff)
    }else{
      yn <- as.numeric(abs(values[hidden2]) < cutoff)
    }
  } else if (sign){
    A <- Ass2
    rownames(A) <- paste0(rownames(Ass2), "_up")
    colnames(A) <- paste0(colnames(Ass2), "_up")
    NA1 <- matrix(0, dim(Ass2)[1], dim(Ass2)[2])
    rownames(NA1) <- paste0(rownames(Ass2), "_up")
    colnames(NA1) <- paste0(colnames(Ass2), "_down")
    NA2 <- matrix(0, dim(Ass2)[1], dim(Ass2)[2])
    rownames(NA2) <- paste0(rownames(Ass2), "_down")
    colnames(NA2) <- paste0(colnames(Ass2), "_up")
    B <- Ass2
    rownames(B) <- paste0(rownames(Ass2), "_down")
    colnames(B) <- paste0(colnames(Ass2), "_down")
    Ass_pm <- rbind(cbind(A, NA1), cbind(NA2, B))
    assign <- apply(Ass_pm,1,function(x) paste(which(x>0)-1,collapse=","))
    terms <- colnames(Ass_pm)
    if(rev){
      yn_p <- as.numeric(abs(values[hidden2]) > cutoff & values[hidden2] > 0)
      yn_m <- as.numeric(abs(values[hidden2]) > cutoff & values[hidden2] < 0)
      yn <- c(yn_p, yn_m)
    }else{
      yn_p <- as.numeric(abs(values[hidden2]) < cutoff & values[hidden2] > 0)
      yn_m <- as.numeric(abs(values[hidden2]) < cutoff & values[hidden2] < 0)
      yn <- c(yn_p, yn_m)
    }
  }
  if(debug){
    print(paste0(dim(assign),dim(terms),dim(yn)))
  }
  write(assign, file=p.assign)
  write.table(terms, file=p.terms, col.names=F, row.names = F, quote=F)
  write.table(yn, file=p.yn, col.names=F, row.names = F, quote=F)
  # running mona (windows/unix)
  if(.Platform$OS.type =="windows"){
    tries <- 0
    while (!file.exists(p.out) & tries<10){
      sys1 <- system(paste( p.mona, "0",
                            p.assign, p.yn, p.terms, p.out, "1",
                            sep = " "), intern = T)
      tries <- tries+1
    }
  } else {
    tries <- 0
    while (!file.exists(p.out) & tries<10){
      sys1 <- system(paste( p.mono, p.mona, "0",
                            p.assign, p.yn, p.terms, p.out, "1",
                            sep = " "), intern = T)
      tries <- tries+1
    }
  }
  
  
  res1 <- read.table(file=p.out, sep="\t", h=F)
  res1 <- res1[order(res1$V2, decreasing=T),]
  colnames(res1) <- c("pathway", "posterior")
  unlink(mona1, recursive = T)
  
  table <- data.frame(res1[, c("pathway", "posterior")], NA, NA, NA)
  colnames(table) <- c("pathway", "score", "score2", "pval", "padj")
  
  res <- list("summary"=table,
              "res"=res1,
              "method"="MONA",
              "opts"=list("data"=data,
                   #       "results"=result,

                          "gmt"=gmt,
                          "minSize"=minSize,
                          "maxSize"=maxSize,
                          "cutoff"=cutoff,
                          "sign"=sign,
                          "rev"=rev))
  
  return(res)
  q(save = "no")
  
}
run_gsea <- function(data, gmt, nperm, minSize, maxSize, debug=FALSE){
  # load libraries
  library(fgsea)
  
  # load data
  if(grepl("\\.RDS", data) | grepl("\\.rds", data)){
    rank <- readRDS(data)
  } else if (grepl("\\.RData", data) | grepl("\\.Rdata", data) | grepl("\\.rdata", data) | grepl("\\.rda", data)){
    # # need to test this further, incomplete (not working)
    # input.env <- environment()
    # load(data, envir = input.env)
    # input.env$data <- ls(envir = input.env) # (not working)
    # rank <- input.env$data
  }
  
  pathways <- gmtPathways(gmt)
  if(debug) head(rank)
  GSEA <- fgsea(pathways,
                rank,
                nperm,
                minSize=minSize,
                maxSize=maxSize)
  
  table <- data.frame(GSEA[, c("pathway", "ES")], NA, GSEA[, c("pval", "padj")])
  colnames(table) <- c("pathway", "score", "score2", "pval", "padj")
  
  res <- list("summary"=table,
              "res"=GSEA,
              "method"="GSEA",
              "opts"=list("data"=data,
                          "results"=GSEA,
                          "gmt"=gmt,
                          "minSize"=minSize,
                          "maxSize"=maxSize,
                          "nperm"=nperm))
  return(res)
}
###################################################################
################### DEBUGGING #####################################
# data <-  "data/mRNA.RDS"
# mygmt <- "data/gmt/c2.cp.kegg.v6.0.symbols.gmt"
# minSize <- 5
# maxSize <- 20
# cutoff  <- 25
# sign <- T
# rev  <- T
# test1 <- run_mona1(data, mygmt, minSize, maxSize, cutoff, sign, rev)
###################################################################