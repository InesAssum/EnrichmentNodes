#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
#' Test script for multiOMICs pathway enrichment
#'
#' @author Kristof Gilicze <kris.gilicze@gmail.com>
#' @param method string // Select method for pathway enrichment
#' @param data string // file path to .RDS input file: named numerical vector, names=gene symbols#
#' @param result string // file path to save results as .RDS 
#' @param gmt string // pathway definition to use: KEGG, GO... default: KEGG
#' @param minSize integer // min size of GSs to be considered, default: 15
#' @param maxSize integer // max size of GSs to be considered, default: 500
#' @param nperm integer // number of permutations, default: 10000
#' @param mygmt string // custom .gmt file // not required
#' 
#' to execute use:
#' Rscript runner.R -f path/to/inputfile -o path/to/output.txt
#' for MONA: Rscript.exe runner.R -m mona -f data/example_data/mRNA.RDS -o test1.RDS -g data/gmt/c2.cp.kegg.v6.0.symbols.gmt --minSize 5 --maxSize 500 -c 25 -s T -r T

# ------------------------------------------------------------------------------
######### library #############
# getwd()
library("optparse")
######### args ################
option_list = list(
  make_option(c("-m", "--method"), type="character", default="gsea", 
              help="Method to use: gsea, mona, mona2d, randomforest (default: gsea)"), 
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name"),
  make_option( "--file2", type="character", default=NULL,
              help="dataset file name of species 2 (mona 2d)"),
  make_option(c("-g", "--gmt"), type="character", default="KEGG",
              help="filepath to gmt pathway definitions (default:KEGG)"),
  make_option("--minSize", default=5, type = "integer",
              help="min size of gene sets considered (default:5)"),
  make_option("--maxSize", type="double", default=500,
              help="max size of gene sets considered (default:500)" ),
  make_option("--nperm", type="double", default=10000,
              help="number of permutations (default: 10000)"),
  make_option(c("-o", "--out"), type="character", default="out.RDS",
              help="output file name [default= out.RDS]"),
  make_option(c("--cutoff", "-c"), type="double", default="25",
              help="MONA: cutoff (default: 25)"),
  make_option("--cutoff2", type="double", default="25",
              help="MONA2d: cutoff species 2 (default: 25)"),
  make_option(c("--sign", "-s"), type="character", default="yes",
              help="MONA: sign (default: yes)"),
  make_option(c("--rev", "-r"), type="character", default="small",
              help="MONA: rev (default: small)"),
  make_option("--debug", type="logical", default="T",
              help="Debug mode (default: F)"),
    make_option("--mygmt", type="character", default="",
              help="Optional: specify your own pathway definition (.gmt).")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

method  <- opt$method
data    <- opt$file
data2   <- opt$file2
out     <- opt$out
gmt     <- opt$gmt
minSize <- opt$minSize
maxSize <- opt$maxSize
nperm   <- opt$nperm
cutoff  <- opt$cutoff
cutoff2 <- opt$cutoff2
rev     <- opt$rev
sign    <- opt$sign
mygmt   <- opt$mygmt
debug   <- opt$debug
################################
######### src ##################
# add methods. Method has to have function named runMETHODNAME.
source("/usr/bin/enrichment_function_library.R")
#getwd()
#print(paste(method, data, gmt, minSize, maxSize, nperm, cutoff, sign, rev, mygmt, out, debug, sep = " "))
################################
if(gmt=="custom"){
	gmt <- mygmt
} else if (gmt=="KEGG"){
  gmt <- "/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"
} else if (gmt=="GO"){
  gmt <- "/data/gmt/c5.bp.v6.0.symbols.gmt"
} else if (gmt=="Reactome"){
  gmt <- "/data/gmt/c2.cp.reactome.v6.0.symbols.gmt"
}
if (method =="gsea"){
  if(debug){
    print(paste(method, data, gmt, minSize, maxSize, nperm, cutoff, sign, rev, mygmt, out, debug, sep = " "))
    getwd()
  }
  result <- run_gsea(data, gmt, nperm, minSize, maxSize, debug)
  saveRDS(result, file = out)
}
if (method =="mona"){
  if(debug){
    print(paste(method, data, gmt, minSize, maxSize, nperm, cutoff, sign, rev, mygmt, out, debug, sep = " "))
    getwd()
  }
  result <- run_mona1(data, gmt, minSize, maxSize, cutoff, sign, rev, debug)
  saveRDS(result, file = out)
}
if (method =="mona2d"){
  if(debug){
    print(paste(method, data, data2, gmt, minSize, maxSize, nperm, cutoff, cutoff2, sign, rev, mygmt, out, debug, sep = " "))
    getwd()
  }
  result <- run_mona2(data, data2, gmt, minSize, maxSize, cutoff, cutoff2, sign, rev, debug)
  saveRDS(result, file = out)
}
if (method =="randomforest"){
  print(paste0(method," is not available. Please try different method."))
}

######### debug ################
# data <- "data/example_data/mRNA.RDS"
# result <- "data/example_data/test1.RDS"
# gmt <- "data/gmt/c2.cp.kegg.v6.0.symbols.gmt"
# minSize <- 5
# maxSize <- 500
# cutoff <- 25
# sign <- T
# rev <- T
################################

