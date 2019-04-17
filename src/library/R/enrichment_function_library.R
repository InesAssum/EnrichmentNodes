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
run_mona1 <- function(data, gmtpath, minSize=5, maxSize=50, cutoff=25, sign="yes", rev="small", debug=FALSE) {
  p.mona <- "/usr/bin/MonaConsoleApp.exe"
  # p.mona <- "/Users/kris.g/maConstruction/software/GKN_plugins/EnrichmentNodesK/src/files/Release/MonaConsoleApp.exe"
  p.mono <- "mono"
  sign <- sign=="yes"
  rev <- rev!="small"
  
  # data reading
  if(debug) print("Reading data...")
  if(grepl("\\.RDS", data) | grepl("\\.rds", data)){
    values <- readRDS(data)
  } else if (grepl("\\.RData", data) | grepl("\\.Rdata", data) | grepl("\\.rdata", data) | grepl("\\.rda", data)){}
  gmt <- gmtPathways(gmtpath)
  if(debug) print("Done.")
  
  # transforming
  if(debug) print("Transforming data...")
  # Assignment matrix
  hidden <- unique(unlist(gmt))
  Ass <- matrix(NA, dimnames = list(hidden, names(gmt)), nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(Ass)[2]){
    Ass[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  hidden2 <- unique(intersect(names(values[!is.na(values)]), hidden))
  if(length(hidden2)==0){stop("GMT file does not match any genes from input. Please select another gmt or check names.")}
  Ass2 <- Ass[hidden2, colnames(Ass)[which(colSums(Ass[hidden2,]) > minSize & colSums(Ass[hidden2,]) < maxSize)]]  # Ass2: Ass, aber nur noch pws groesser und kleiner als min max
  if(debug) print("Done.")

  # temp creation  
  if(debug) print("Creating temp files...")
  mona1 <- tempdir()
  dir.create(mona1, showWarnings = F)
  if(debug) print(mona1)
  p.out <- paste0(mona1, "/", "output.txt")
  p.assign <- paste0(mona1, "/", "assignmentMatrix.txt")
  p.terms <- paste0(mona1, "/", "terms.txt")
  p.yn <- paste0(mona1, "/", "values.txt")
  if(debug) print("Done.")
  
  # transform and filter
  if(debug) print("Filtering...")
  if(!sign){      #  Richtung der regulierung, negative = runter reguliert
    assign <- apply(Ass2,1,function(x) paste(which(x>0)-1,collapse=","))
    terms <- colnames(Ass2)
    if(rev){
      yn <- as.numeric(abs(values[hidden2]) > cutoff)
    }else{
      yn <- as.numeric(abs(values[hidden2]) < cutoff)
    }
  } else if (sign){
    # Assignment matrix distinguishing values >/< 0
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
  # storing inbetween files
  write(assign, file=p.assign)
  write.table(terms, file=p.terms, col.names=F, row.names = F, quote=F)
  write.table(yn, file=p.yn, col.names=F, row.names = F, quote=F)
  # running mona (windows/unix)
  if(debug) print("Done.")
  if(debug) print("Running mona...")
  if(.Platform$OS.type =="windows"){
    tries <- 0
    while (!file.exists(p.out) && tries<3){
      sys1 <- system(paste( p.mona, "0",
                            p.assign, p.yn, p.terms, p.out, "1",
                            sep = " "), intern = T)
      tries <- tries+1
    }
  } else {
    tries <- 0
    while (!file.exists(p.out) && tries<3){
      sys1 <- system(paste( p.mono, p.mona, "0",
                            p.assign, p.yn, p.terms, p.out, "1",
                            sep = " "), intern = T)
      tries <- tries+1
    }
  }
  if(tries==3) stop("Failed to execute MonaConsoleApp.exe.")
  if(debug) print("Done.")
  if(debug) print("Creating summary...")
  
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
  if(debug) print("Done.")
  return(res)
  q(save = "no")
  
}
#' Mona 2D
run_mona2 <- function(data1, data2, gmtpath, minSize=5, maxSize=50, cutoff1=25, cutoff2=25, sign="yes", rev="small", debug=FALSE) {
  p.mona <- "/usr/bin/MonaConsoleApp.exe"
  # p.mona <- "/Users/kris.g/maConstruction/software/GKN_plugins/EnrichmentNodesK/src/files/Release/MonaConsoleApp.exe"
  p.mono <- "mono"
  sign <- sign=="yes"
  rev <- rev!="small"
  values1=NULL
  values2=NULL
  
  # data reading
  if(debug) print("Reading data...")
  try(if(grepl("\\.RDS", data1) | grepl("\\.rds", data1)){
    values1 <- readRDS(data1)
    values2 <- readRDS(data2)
  })
  try(if (typeof(data1)=="double"){
    values1 <- data1
    values2 <- data2
  })
  if(is.null(values1) || is.null(values2)) stop("Can't deal with this input data")
  
  gmt <- gmtPathways(gmtpath)
  if(debug) print("Done.")
  
  # transforming
  if(debug) print("Transforming data...")
  # Assignment matrix
  hidden <- unique(unlist(gmt))
  Ass <- matrix(NA, dimnames = list(hidden, names(gmt)), nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(Ass)[2]){
    Ass[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  hidden2 <- unique(intersect(names(values1[!is.na(values1)]), hidden))
  if(length(hidden2)==0){stop("GMT file does not match any genes from species1. Please select another gmt or check names.")}
  Ass2 <- Ass[hidden2, colnames(Ass)[which(colSums(Ass[hidden2,]) > minSize & colSums(Ass[hidden2,]) < maxSize)]]  # Ass2: Ass mit pws groesser und kleiner als min max
  
  if(debug) print("Done.")
  
  # temp creation  
  if(debug) print("Creating temp files...")
  mona1 <- tempdir()
  dir.create(mona1, showWarnings = F)
  if(debug) print(mona1)
  p.out <- paste0(mona1, "/", "output.txt")
  p.assign <- paste0(mona1, "/", "assignmentMatrix.txt")
  p.terms <- paste0(mona1, "/", "terms.txt")
  p.species1 <- paste0(mona1, "/", "species1.txt")
  p.species2 <- paste0(mona1, "/", "species2.txt")
  p.species2_miss <- paste0(mona1, "/", "species2_miss.txt")
  if(debug) print("Done.")
  
  # transform and filter
  if(debug) print("Filtering...")
  if(!sign){      #  sign = Richtung der regulierung, negative = runter reguliert...
    assign <- apply(Ass2,1,function(x) paste(which(x>0)-1,collapse=","))
    terms <- colnames(Ass2)
    species2_miss <- as.numeric(!(hidden2 %in% names(values2)))
    if(rev){
      species1 <- as.numeric(abs(values1[hidden2]) > cutoff1)
      species2 <- as.numeric(abs(values2[hidden2]) > cutoff2)
      species2[is.na(species2)] <- 0
    }else{
      species1 <- as.numeric(abs(values1[hidden2]) < cutoff1)
      species2 <- as.numeric(abs(values2[hidden2]) < cutoff2)
      species2[is.na(species2)] <- 0
    }
  } else if (sign){
    # Assignment matrix distinguishing values >/< 0
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
    species2_miss <- c(as.numeric(!(hidden2 %in% names(values2))),as.numeric(!(hidden2 %in% names(values2))))
    
    
    if(rev){
      species1_p <- as.numeric(abs(values1[hidden2]) > cutoff1 & values1[hidden2] > 0)
      species1_m <- as.numeric(abs(values1[hidden2]) > cutoff1 & values1[hidden2] < 0)
      species1 <- c(species1_p, species1_m)
      
      species2_p <- as.numeric(abs(values2[hidden2]) > cutoff2 & values2[hidden2] > 0)
      species2_m <- as.numeric(abs(values2[hidden2]) > cutoff2 & values2[hidden2] < 0)
      species2 <- c(species2_p, species2_m)
      species2[is.na(species2)] <- 0
    }else{
      species1_p <- as.numeric(abs(values1[hidden2]) < cutoff1 & values1[hidden2] > 0)
      species1_m <- as.numeric(abs(values1[hidden2]) < cutoff1 & values1[hidden2] < 0)
      species1 <- c(species1_p, species1_m)
      
      species2_p <- as.numeric(abs(values2[hidden2]) < cutoff2 & values2[hidden2] > 0)
      species2_m <- as.numeric(abs(values2[hidden2]) < cutoff2 & values2[hidden2] < 0)
      species2 <- c(species2_p, species2_m)
      species2[is.na(species2)] <- 0
    }
  }

  # storing inbetween files
  write(assign, file=p.assign)
  write.table(terms, file=p.terms, col.names=F, row.names = F, quote=F)
  write.table(species1, file=p.species1, col.names=F, row.names = F, quote=F)
  write.table(species2, file=p.species2, col.names=F, row.names = F, quote=F)
  write.table(species2_miss, file=p.species2_miss, col.names=F, row.names = F, quote=F)
  # running mona (windows/unix)
  if(debug) print("Done.")
  if(debug) print("Running mona...")

    tries <- 0
    while (!file.exists(p.out) && tries<3){
      sys1 <- system(paste( p.mono, p.mona, "7",
                            p.assign, p.species1, p.terms, p.out, "1", p.species2, p.species2_miss,
                            sep = " "), intern = T)
      tries <- tries+1
    }
  
  if(tries==3) stop("Failed to execute MonaConsoleApp.exe.")
  if(debug) print("Done.")
  if(debug) print("Creating summary...")
  
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
                          "cutoff1"=cutoff1,
                          "cutoff2"=cutoff2,
                          "sign"=sign,
                          "rev"=rev))
  if(debug) print("Done.")
  return(res)
  q(save = "no")
  
}
run_gsea <- function(data, gmtpath, nperm=1000, minSize=5, maxSize=50, debug=FALSE){
  # load data
  if(debug) print("Reading data...")
  if(grepl("\\.RDS", data) | grepl("\\.rds", data)){
    rank <- readRDS(data)
  } else if (grepl("\\.RData", data) | grepl("\\.Rdata", data) | grepl("\\.rdata", data) | grepl("\\.rda", data)){
    # # need to test this further, incomplete (not working)
    # input.env <- environment()
    # load(data, envir = input.env)
    # input.env$data <- ls(envir = input.env) # (not working)
    # rank <- input.env$data
  }
  
  pathways <- gmtPathways(gmtpath)
  if(length(intersect(unique(unlist(pathways)), names(rank)))==0){stop("GMT file does not match any genes from input. Please select another gmt or check names.")}
  if(debug) print("Done.")
  if(debug) print("Running fgsea...")
  GSEA <- fgsea(pathways,
                rank,
                nperm,
                minSize=minSize,
                maxSize=maxSize)
  if(debug) print("Done.")
  if(debug) print("Creating result...")
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
  if(debug) print("Done.")
  return(res)
}
###################################################################
################### DEBUGGING #####################################
  # data <-  "/Users/kris.g/maConstruction/scripts/data/example_data/mRNA.RDS"
  # gmtpath <- "/Users/kris.g/maConstruction/scripts/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"
  # gmtpath <- "/Users/kris.g/Documents/gmts/drekegg.gmt"
  # test1 <- run_mona1(data, gmtpath, sign="no",rev = "big")
  # test2 <- run_gsea(data, gmtpath)
  # values2 <- values1[1:200]
  # saveRDS(values2, "/Users/kris.g/maConstruction/scripts/data/example_data/prot.RDS")
  # data1 <- "/Users/kris.g/maConstruction/scripts/data/example_data/mRNA.RDS"
  # data2 <- "/Users/kris.g/maConstruction/scripts/data/example_data/prot.RDS"
  # test3 <- run_mona2(data1, data2, gmtpath, sign="no",rev = "big")
  ###################################################################
