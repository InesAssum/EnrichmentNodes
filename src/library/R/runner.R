#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
#' Scripts for multiOMICs pathway enrichment
#'
#' @param method string // Select method for pathway enrichment
#' @param data string // file path to .RDS input file: named numerical vector, names=gene symbols#
#' @param result string // file path to save results as .RDS 
#' @param gmt string // pathway definition to use: KEGG, GO... default: KEGG
#' @param minSize integer // min size of GSs to be considered, default: 15
#' @param maxSize integer // max size of GSs to be considered, default: 500
#' @param type string // type of result information, e.g. pvalue, score or significant
#' @param cutoff numeric // cutoff value for significance
#' @param sign string // taking sign into account yes or no
#' @param nperm integer // Deprecated: number of permutations, default: 0
#' @param twomiss string // running MONA2Dcoop with two missings, default: no
#' @param mygmt string // custom .gmt file // not required
#' 
#' @author Ines Assum (previous work from Kristof Gilicze included)
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
  make_option(c("-m", "--method"), type="character", default="GSEA", 
              help="Method to use: GSEA, MGSA, singleMONA, MONA2Dcoop, MONA2Dinh, MONA3Dcoop, randomforest (default: gsea)"), 
  make_option(c("-f", "--file"), type="character", default="/data/example_data/mRNA.xlsx",
              help="dataset file name (first/single species)"),
  make_option("--file1", type="character", default="/data/example_data/mRNA.xlsx",
              help="dataset file name (first/single species)"),
  make_option( "--file2", type="character", default="/data/example_data/prot.xlsx",
              help="dataset file name of species 2 (MONA2Dcoop/MONA2Dinh)"),
  make_option( "--file3", type="character", default="/data/example_data/meta.xlsx",
               help="dataset file name of species 3 (MONA3Dcoop)"),
  make_option(c("-g", "--gmt"), type="character", default="KEGG",
              help="filepath to gmt pathway definitions (default:KEGG)"),
  make_option("--minSize", default=5, type = "integer",
              help="min size of gene sets considered (default:5)"),
  make_option("--maxSize", type="double", default=500,
              help="max size of gene sets considered (default:500)" ),
  make_option(c("--type", "-t"), type="character", default="pvalue",
              help="GSEA/MGSA/MONA type of ranking, pvalue, score or significant (default: pvalue)"),
  make_option(c("--cutoff", "-c"), type="double", default="0.05",
              help="MGSA/MONA: Significance: cutoff (default: 0.05)"),
  make_option(c("--sign", "-s"), type="character", default="no",
              help="GSEA/MGSA/MONA: take direction of effect into account? (default: no)"),
  make_option("--nperm", type="double", default=0,
              help="Deprecated: GSEA number of permutations (default: NULL)"),
  make_option(c("--twomiss"), type="character", default="no",
              help="Consider MONA2D with missings for both omics (default: no)"),
  # make_option("--cutoff2", type="double", default="25",
  #             help="MONA2d: cutoff species 2 (default: 25)"),
  make_option(c("-o", "--outtable"), type="character", default="out.txt",
              help="output file name [default= out.txt]"),
  make_option(c("--outres"), type="character", default="out.RDS",
              help="output file name [default= out.RDS]"),
  make_option("--debug", type="logical", default=F,
              help="Debug mode (default: F)"),
  make_option("--mygmt", type="character", default=NULL,
              help="Optional: specify your own pathway definition (.gmt).")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

method  <- opt$method
data    <- opt$file
data1   <- opt$file1
data2   <- opt$file2
data3   <- opt$file3
outres  <- opt$outres
outtxt  <- opt$outtable
gmt     <- opt$gmt
minSize <- opt$minSize
maxSize <- opt$maxSize
type    <- opt$type
cutoff  <- opt$cutoff
sign    <- opt$sign
twomiss <- opt$twomiss
nperm   <- opt$nperm
#cutoff2 <- opt$cutoff2
mygmt   <- opt$mygmt
debug   <- opt$debug
################################
######### src ##################
# add methods. Method has to have function named runMETHODNAME.
source("/usr/bin/enrichment_function_library.R")
#getwd()
#print(paste(method, data, gmt, minSize, maxSize, nperm, cutoff, sign, rev, mygmt, outres, outtxt, debug, sep = " "))
################################
if(gmt=="custom"){
	gmt <- mygmt
} else if (gmt=="KEGG"){
  gmt <- "/data/gmt/c2.cp.kegg.v7.5.1.symbols.gmt"
} else if (gmt=="GO:BP"){
  gmt <- "/data/gmt/c5.go.bp.v7.5.1.symbols.gmt"
} else if (gmt=="Reactome"){
  gmt <- "/data/gmt/c2.cp.reactome.v7.5.1.symbols.gmt"
} else if (gmt=="CPDB"){
  gmt <- "/data/gmt/CPDB_pathways_genes_metabolites.gmt"
}
if (method =="GSEA"){
  if(debug){
    print(paste(method, data, gmt, minSize, maxSize, type, cutoff, sign, nperm, mygmt, outres, outtxt, debug, sep = " "))
    getwd()
  }
  result <- run_gsea(data, gmt, minSize, maxSize, type, cutoff, sign, nperm, debug)
  saveRDS(result, file = outres)
  write.table(result$summary, file = outtxt,
              row.names = F, col.names = T, quote = F, sep = "\t")
}
if (method =="MGSA"){
  if(debug){
    print(paste(method, data, gmt, minSize, maxSize, type, cutoff, sign, mygmt, outres, outtxt, debug, sep = " "))
    getwd()
  }
  result <- run_mgsa(data, gmt, minSize, maxSize, type, cutoff, sign, debug)
  saveRDS(result, file = outres)
  write.table(result$summary, file = outtxt,
              row.names = F, col.names = T, quote = F, sep = "\t")
}
if (method =="singleMONA"){
  if(debug){
    print(paste(method, data, gmt, minSize, maxSize, type, cutoff, sign, mygmt, outres, outtxt, debug, sep = " "))
    getwd()
  }
  result <- run_single_mona(data, gmt, minSize, maxSize, type, cutoff, sign, debug)
  saveRDS(result, file = outres)
  write.table(result$summary, file = outtxt,
              row.names = F, col.names = T, quote = F, sep = "\t")
}
if (method =="MONA2Dcoop"){
  if(debug){
    print(paste(method, data1, data2, gmt, minSize, maxSize, type, cutoff, sign, mygmt, outres, outtxt, debug, sep = " "))
    getwd()
  }
  if(twomiss=="yes"){
    result <- run_mona2Dcoop2(data1, data2, gmt, minSize, maxSize, type, cutoff, sign, debug)
  }else{
    result <- run_mona2Dcoop(data1, data2, gmt, minSize, maxSize, type, cutoff, sign, debug)
  }
  saveRDS(result, file = outres)
  write.table(result$summary, file = outtxt,
              row.names = F, col.names = T, quote = F, sep = "\t")
}
if (method =="MONA2Dinh"){
  print(paste(method, data1, data2, gmt, minSize, maxSize, type, cutoff, sign, mygmt, outres, outtxt, debug, sep = " "))
  getwd()
  print(paste0(method," is not available yet. Please try different method."))
  # result <- run_mona2Dinh(data1, data2, gmt, minSize, maxSize, type, cutoff, sign, debug)
  # saveRDS(result, file = outres)
  # write.table(result$summary, file = outtxt,
  #             row.names = F, col.names = T, quote = F, sep = "\t")
}
if (method =="MONA3Dcoop"){
  if(debug){
    print(paste(method, data1, data2, data3, gmt, minSize, maxSize, type, cutoff, sign, mygmt, outres, outtxt, debug, sep = " "))
    getwd()
  }
    result <- run_mona3Dcoop(data1, data2, data3, gmt, minSize, maxSize, type, cutoff, sign, debug)
  saveRDS(result, file = outres)
  write.table(result$summary, file = outtxt,
              row.names = F, col.names = T, quote = F, sep = "\t")
}
if (method =="randomforest"){
  print(paste0(method," is not available yet. Please try different method."))
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

