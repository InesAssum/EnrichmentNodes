#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ------------------------------------------------------------------------------
#' Test script to run MONA from commandline using generic KNIME nodes
#'
#' @author Ines Assum
#' @param data string // file path to .RDS input file: named numerical vector, names=gene symbols#
#' @param result string // file path to save results as .RDS 
#' @param gmt string // Gene set definition to be used (KEGG/GO/custom), default: KEGG
#' @param minSize integer // min size of GSs to be considered, default: 15
#' @param maxSize integer // max size of GSs to be considered, default: 500
#' @param cutoff numeric // cutoff for significance threshold, default: 0.05
#' @param sign bool // integrate direction of effect? default: 0
#' @param rev bool // pvalues (small values sign) vs. counts (big values significant), default: 0
#' @param mygmt string // custom .gmt file // not necessary
#' 
#' to execute use:
#' run_mona1.R 'data.path' 'out.path' 'KEGG' '15' '500' '0.05' '0' '0' 'customgmtfile'
# ------------------------------------------------------------------------------

data <- args[1]
data
result <- args[2]
result
gmt2 <- args[3]
gmt2
minSize <- as.integer(args[4])
minSize
maxSize <- as.integer(args[5])
maxSize
cutoff <- as.numeric(args[6])
cutoff
sign <- as.logical(args[7]=="yes")
args[7]
sign
rev <- !as.logical(args[8]=="small")
args[8]
rev
mygmt <- args[9]
mygmt

#./usr/bin/run_mona_single.R '/data/example_data/mRNA.RDS' '/data/example_data/test1.RDS' 'KEGG' '15' '500' '25' 'yes' 'large' '/data/gmt/c2.cp.kegg.v6.0.symbols.gmt'

# data <- "/data/example_data/mRNA.RDS"
# result <- "/data/example_data/test1.RDS"
# gmt2 <- "custom"
# minSize <- 5
# maxSize <- 500
# cutoff <- 25
# sign <- T
# rev <- T
# mygmt <- "/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"

# data <- "/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/software/GKN_plugins/EnrichmentNodes/src/data/mRNA.RDS"
# result <- "/Users/ines/Downloads/test_local.RDS"
# gmt2 <- "custom"
# minSize <- 5
# maxSize <- 500
# cutoff <- 25
# sign <- "true"
# rev <- "true"
# mygmt <- "/Users/ines/Documents/ICB/PhD/projects/KNIME_multiOMICs_enrichment/software/GKN_plugins/EnrichmentNodes/src/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"

# data <- "/home/icb/ines.assum/projects/KNIME_multiOMICs_enrichment/software/GKN_plugins/EnrichmentNodes/src/data/mRNA.RDS"
# result <- "/home/icb/ines.assum/rstudio/test_server.RDS"
# gmt2 <- "custom"
# minSize <- 5
# maxSize <- 500
# cutoff <- 25
# sign <- "true"
# rev <- "true"
# mygmt <- "/home/icb/ines.assum/projects/KNIME_multiOMICs_enrichment/software/GKN_plugins/EnrichmentNodes/src/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"


if(gmt2=="KEGG"){gmt2 <- "/data/gmt/c2.cp.kegg.v6.0.symbols.gmt"}
if(gmt2=="GO"){gmt2 <- "/data/gmt/c5.bp.v6.0.symbols.gmt"}
if(gmt2=="custom"){gmt2 <- mygmt}
gmt2

p.mono <- "mono"
p.mona <- "/usr/bin/MonaConsoleApp.exe"
# p.mona <- "/Users/ines/Documents/ICB/PhD/projects/MONA/software/MonaConsoleApp_3.1/MonaConsoleApp/bin/Debug/MonaConsoleApp.exe"
# p.mona <- "/Users/ines/Documents/ICB/PhD/projects/MONA/software/MonaConsoleApp_Andi_mod/MonaConsoleApp/bin/Debug/MonaConsoleApp.exe"
# p.mono <- "cd /home/icb/ines.assum/ext_tools/mono/bin \n ./mono"
# p.mona <- "/home/icb/ines.assum/projects/MONA/software/MonaConsoleApp_Andi_mod/MonaConsoleApp/bin/Debug/MonaConsoleApp.exe"
# p.mona <- "/home/icb/ines.assum/projects/KNIME_multiOMICs_enrichment/software/GKN_plugins/EnrichmentNodes/docker/enrich/files/MonaConsoleApp.exe"



library(fgsea)

gmt <- gmtPathways(gmt2)
hidden <- unique(unlist(gmt))
Ass <- matrix(NA, dimnames = list(hidden, names(gmt)), nrow = length(hidden), ncol = length(gmt))
for (i in 1:dim(Ass)[2]){
  Ass[,i] <- as.numeric(hidden %in% gmt[[i]])
}

if(grepl("\\.RDS", data) | grepl("\\.rds", data)){
  values <- readRDS(data)
} else if (grepl("\\.RData", data) | grepl("\\.Rdata", data) | grepl("\\.rdata", data) | grepl("\\.rda", data)){
  # # need to test this further, incomplete (not working)
  # input.env <- environment()
  # load(data, envir = input.env)
  # input.env$data <- ls(envir = input.env) # (not working)
  # rank <- input.env$data
}
head(values)
hidden2 <- unique(intersect(names(values[!is.na(values)]),
                            hidden))

length(hidden2)
Ass2 <- Ass[hidden2, colnames(Ass)[which(colSums(Ass[hidden2,]) > minSize &
                                           colSums(Ass[hidden2,]) < maxSize)]]

mona1 <- tempdir()
mona1
p.out <- paste0(mona1, "/", "output.txt")
p.assign <- paste0(mona1, "/", "assignmentMatrix.txt")
p.terms <- paste0(mona1, "/", "terms.txt")
p.yn <- paste0(mona1, "/", "values.txt")

str(sign)

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

write(assign, file=p.assign)
write.table(terms, file=p.terms, col.names=F, row.names = F, quote=F)
write.table(yn, file=p.yn, col.names=F, row.names = F, quote=F)


tries <- 0
while (!file.exists(p.out) & tries<10){
  sys1 <- system(paste(p.mono, p.mona, "0",
                       p.assign, p.yn, p.terms, p.out, "1",
                       sep = " "), intern = T)
  tries <- tries+1
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
                        "results"=result,
                        "gmt"=gmt,
                        "minSize"=minSize,
                        "maxSize"=maxSize,
                        "cutoff"=cutoff,
                        "sign"=sign,
                        "rev"=rev))

saveRDS(res, file=result)
q(save = "no")







