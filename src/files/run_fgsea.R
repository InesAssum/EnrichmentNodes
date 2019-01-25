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
minSize <- as.integer(args[4])
maxSize <- as.integer(args[5])
nperm <- as.integer(args[6])
mygmt <- args[7]
#mygmt

# data <- "/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/software/GKN_plugins/EnrichmentNodes/src/data/mRNA.RDS"
# result <- "/Users/ines/Downloads/test_local.RDS"
# gmt <- "custom"
# minSize <- 5
# maxSize <- 500
# nperm <- 10000
# mygmt <- "/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/software/GKN_plugins/EnrichmentNodes/src/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"

# data <- "/data/example_data/mRNA.RDS"
# result <- "/data/example_data/test1.RDS"
# gmt <- "custom"
# minSize <- 5
# maxSize <- 500
# nperm <- 10000
# mygmt <- "/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"


if(gmt=="KEGG"){gmt <- "/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"}
if(gmt=="GO"){gmt <- "/data/gmt/c5.bp.v6.0.symbols.gmt"}
if(gmt=="custom"){gmt <- mygmt}

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
head(rank)
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
                        "results"=result,
                        "gmt"=gmt,
                        "minSize"=minSize,
                        "maxSize"=maxSize,
                        "nperm"=nperm))

saveRDS(res, file=result)
q(save = "no")

  