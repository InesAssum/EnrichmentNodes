#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ------------------------------------------------------------------------------
#' Test script to call GSEA analysis from commandline using generic KNIME nodes
#'
#' @author Ines Assum
#' @param data string // file path to .RDS input file: named numerical vector, names=gene symbols#
#' @param result string // file path to save results as .RDS 
#' @param gmt string // Gene set definition to be used (KEGG/GO/custom), default: KEGG
#' @param minSize integer // min size of GSs to be considered, default: 15
#' @param maxSize integer // max size of GSs to be considered, default: 500
#' @param nperm integer // number of permutations, default: 10000
#' @param mygmt string // custom .gmt file // not necessary
#' 
#' to execute use:
#' #R CMD BATCH '--args data="/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/local_ines/example_data/mRNA.RDS" result="/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/local_ines/example_data/test_result.RDS" gmt="/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/local_ines/example_data/gmt/c2.cp.kegg.v6.0.symbols.gmt" minSize=5 maxSize=500 nperm=10000' GKN_run_fgsea.R &
#' run_fgsea.R 'data.path' 'out.path' 'KEGG' '15' '500' '10000' 'customgmtfile'
# ------------------------------------------------------------------------------

#' Test script to call GSEA analysis from commandline using generic KNIME nodes
#' 
#' 
# # parse arguments - deprecated - use positional arguments
# args=(commandArgs(TRUE))
# eval(parse(text=args))

# data: path to .RDS data file
# result: path to save results as .RDS file

data <- args[1]
result <- args[2]
gmt <- args[3]
minSize <- args[4]
maxSize <- args[5]
nperm <- args[6]
mygmt <- args[7]

if(gmt=="KEGG"){gmt <- "/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"}
if(gmt=="GO"){gmt <- "/data/gmt/c5.bp.v6.0.symbols.gmt"}
if(gmt=="custom"){gmt <- mygmt}


# load libraries
library(fgsea)

# load data
if(grep("\\.RDS | \\.rds", data)){
  rank <- readRDS(data)
} else if (grep("\\.RData | \\.Rdata | \\.rdata", data)){
  # # need to test this further, incomplete (not working)
  # input.env <- environment()
  # load(data, envir = input.env)
  # input.env$data <- ls(envir = input.env) # (not working)
  # rank <- input.env$data
}

pathways <- gmtPathways(gmt)
GSEA <- fgsea(pathways = pathways,
              stats = rank,
              minSize=minSize,
              maxSize=maxSize,
              nperm=nperm)
table[, c("Pathway", "score", "score2", "pvalue", "padjust")] <-
  c(GSEA[, c("Pathway", "ES")], NA, GSEA[, c("pvalue", "padjust")])

res <- list("summary"=table,
            "res"=GSEA,
            "method"="GSEA",
            "opts"=list("data"=data,
                        "results"=result,
                        "gmt"=gmt,
                        "minSize"=minSize,
                        "maxSize"=maxSize,
                        "nperm"=nperm))

saveRDS(res, file=result)
q(save = "no")

  