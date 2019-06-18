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
    if(!rev){
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
    if(!rev){
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
run_mona2 <- function(data1, data2, gmtpath, minSize=5, maxSize=50, cutoff1=25, cutoff2=25, sign="yes", rev="small", debug=FALSE) {p.mona <- "/usr/bin/MonaConsoleApp.exe"
p.mona <- "/usr/bin/MonaConsoleApp.exe"
# p.mona <- "/Users/kris.g/Documents/old_mona_stuff/MonaConsoleApp.exe"
p.mono <- "mono"
sign <- sign=="yes"
rev <- rev!="small"

# data reading
if(debug) print("Reading data...")
try(if(typeof(data1)=="double"){
  values1 <- data1; values2 <- data2
} 
else if(grepl("\\.RDS", data1) | grepl("\\.rds", data1)){
  values1 <- readRDS(data1); values2 <- readRDS(data2)
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
hidden2 <- unique(intersect(c(names(values1[!is.na(values1)]),names(values2[!is.na(values2)])), hidden))
if(length(hidden2)==0){stop("GMT file does not match any genes from species1. Please select another gmt or check names.")}
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
p.species1 <- paste0(mona1, "/", "species1.txt")
p.species2 <- paste0(mona1, "/", "species2.txt")
p.species1_miss <- paste0(mona1, "/", "species1_miss.txt") 
p.species2_miss <- paste0(mona1, "/", "species2_miss.txt")
if(debug) print("Done.")

# transform and filter
if(debug) print("Filtering...")
if(!sign){      #  sign = Richtung der regulierung, negative = runter reguliert...
  assign <- apply(Ass2,1,function(x) paste(which(x>0)-1,collapse=","))
  terms <- colnames(Ass2)
  
  species1_miss <- as.numeric(!(hidden2 %in% names(values1)))
  species2_miss <- as.numeric(!(hidden2 %in% names(values2)))
  
  species1 <- as.numeric(abs(values1[hidden2]) > cutoff1)
  species1[is.na(species1)] <- 0
  species2 <- as.numeric(abs(values2[hidden2]) > cutoff2)
  species2[is.na(species2)] <- 0
} else if (sign){
  # stop("MONA2D extended is not supported with sign = TRUE.")
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
  
  species1_miss <- c(as.numeric(!(hidden2 %in% names(values1))),as.numeric(!(hidden2 %in% names(values1))))
  species2_miss <- c(as.numeric(!(hidden2 %in% names(values2))),as.numeric(!(hidden2 %in% names(values2))))
  species1_p <- as.numeric(abs(values1[hidden2]) > cutoff1 & values1[hidden2] > 0)
  species1_m <- as.numeric(abs(values1[hidden2]) > cutoff1 & values1[hidden2] < 0)
  species1 <- c(species1_p, species1_m)
  species1[is.na(species1)] <- 0
  
  species2_p <- as.numeric(abs(values2[hidden2]) > cutoff2 & values2[hidden2] > 0)
  species2_m <- as.numeric(abs(values2[hidden2]) > cutoff2 & values2[hidden2] < 0)
  species2 <- c(species2_p, species2_m)
  species2[is.na(species2)] <- 0
}
if(rev){
  species1 <- as.numeric(species1==0)
  species2 <- as.numeric(species2==0)
}
# storing inbetween files
write(assign, file=p.assign)
write.table(terms, file=p.terms, col.names=F, row.names = F, quote=F)
write.table(species1, file=p.species1, col.names=F, row.names = F, quote=F)
write.table(species2, file=p.species2, col.names=F, row.names = F, quote=F)
write.table(species1_miss, file=p.species1_miss, col.names=F, row.names = F, quote=F)
write.table(species2_miss, file=p.species2_miss, col.names=F, row.names = F, quote=F)
# running mona (windows/unix)
if(debug) print("Done.")
if(debug) print("Running mona...")

tries <- 0
while (!file.exists(p.out) && tries<3){
  sys1 <- system(paste( p.mono, p.mona, "9",
                        p.assign, p.species1, p.terms, p.out, "1", p.species2, p.species1_miss, p.species2_miss,
                        sep = " "), intern = T)
  tries <- tries+1
}

if(tries==3) stop("Failed to execute MonaConsoleApp.exe.")
if(debug) print("Done.")
if(debug) print("Creating summary...")

res1 <- read.table(file=p.out, sep="\t", h=F, stringsAsFactors = F)
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
#' FGSEA
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
###################################################################
################# Simulation and Evaluation #######################
#' Z Score Simulation
simulate_zscores <- function(gmtpath, N=5,rho=0.5, coverage=0.25, minsize=10, include_direction=T, dynamic_SD=F, sdsig=5){
  print("Initialising...")
  ## Variables
  # Correlation coefficient rho
  # Number of simulated data sets N
  # Min size of pathway minsize
  # Error rates
  fpr <- c(1e-4,c((1:5)/10))
  fnr <- fpr
  
  ## Read:
  # GMT
  if(gmtpath=="mac") gmtpath <- "/Users/kris.g/maConstruction/scripts/SimEval/c2.cp.kegg.v6.0.symbols.gmt"   ### <- for debugging
  gmt <- gmtPathways(gmtpath)
  genes.kegg <- unique(unlist(gmt))
  pathways <- names(gmt)
  ass.kegg <- matrix(NA, dimnames = list(genes.kegg, names(gmt)),nrow = length(genes.kegg), ncol = length(gmt))
  for (i in 1:ncol(ass.kegg)){ ass.kegg[,i] <- as.numeric(genes.kegg %in% gmt[[i]])}     # genes x pathways
  
  # Trans and Prot with cov 25%
  genes.t <- sample(genes.kegg, round(length(genes.kegg)*0.8))
  genes.p <- c(setdiff(genes.kegg, genes.t), sample(genes.t, round(coverage*length(genes.t))))
  # length(intersect(genes.t,genes.p))/length(genes.t)
  genes.tp <- unique(c(genes.t, genes.p))
  
  
  # Assignment matrices
  ass.trans <- ass.kegg[genes.t, colnames(ass.kegg)[colSums(ass.kegg[genes.t,])>=minsize]]
  ass.prot <- ass.kegg[genes.p, colnames(ass.kegg)[colSums(ass.kegg[genes.p,])>=minsize]]
  pw.trans <- colnames(ass.trans)
  pw.prot <- colnames(ass.prot)
  
  ## Init
  # Z scores
  sim.data <- list()
  
  # Options
  with_sign = include_direction
  
  
  print("Starting Simulation...")
  print(paste("Estimated time:",.ntupel(fpr)*N*0.04,"seconds."))
  
  
  for(k in (1:N)){
    print(paste("Run",k,"of",N))
    for(alpha in fpr){
      for(beta in fnr){
        if(alpha+beta>1) break
        
        # Normal distribution of background and signal (+/-)
        if(dynamic_SD){
          sd.bg <- 1+100*beta
          sd.sig <- (1-alpha)*10
        } else {
          sd.bg <- 1  
          sd.sig <- sdsig
        }
        cutoff.a <- qnorm(1-alpha/2, sd=sd.bg)
        cutoff.b <- qnorm(beta, sd=sd.sig)
        mean.bg <- 0
        mean.sig <- cutoff.a - cutoff.b
        
        # Covariance matrix Sigma for bivariate distributen
        sigma.bg <- matrix(c(sd.bg^2, sd.bg^2*rho, sd.bg^2*rho, sd.bg^2), 2)
        sigma.tp <- matrix(c(sd.sig^2, sd.sig^2*rho, sd.sig^2*rho, sd.sig^2), 2)
        sigma.t <- matrix(c(sd.sig^2, sd.sig*sd.bg*rho, sd.sig*sd.bg*rho, sd.bg^2), 2)
        sigma.p <- matrix(c(sd.bg^2, sd.bg*sd.sig*rho, sd.bg*sd.sig*rho, sd.sig^2), 2)
        
        # Sample active terms and their sign (3 per transcriptome, 3 per proteome
        act.terms <- data.frame(terms.trans=as.character(sample(pw.trans, 3)),
                                signs.trans=NA,
                                terms.prot=as.character(sample(pw.prot, 3)),
                                signs.prot=NA,
                                stringsAsFactors = F)
        active.terms <- unique(c(act.terms$terms.trans,
                                 act.terms$terms.prot))
        active.terms <- data.frame(term=active.terms,
                                   sign=sample(c(-1,1), length(active.terms),replace=T),
                                   row.names = active.terms,
                                   stringsAsFactors = F)
        
        act.terms$signs.trans <- active.terms[act.terms$terms.trans, "sign"]
        act.terms$signs.prot <- active.terms[act.terms$terms.prot, "sign"]
        
        # Active genes    
        act.genes.trans <- data.frame(gene = rownames(ass.trans)[rowSums(ass.trans[genes.t, act.terms$terms.trans])>0],
                                      sign = NA,
                                      stringsAsFactors = F)
        rownames(act.genes.trans) <- act.genes.trans$gene
        act.genes.prot <- data.frame(gene = rownames(ass.prot)[rowSums(ass.prot[genes.p, act.terms$terms.prot])>0],
                                     sign = NA,
                                     stringsAsFactors = F)
        rownames(act.genes.prot) <- act.genes.prot$gene
        active.genes <- data.frame(gene = as.character(unique(c(act.genes.trans$gene, act.genes.prot$gene))),
                                   sign = NA,
                                   stringsAsFactors = F)
        rownames(active.genes) <- active.genes$gene
        # sample sign, in case it belongs to multiple active pathways
        for (l in 1:length(active.genes$gene)){
          active.genes$sign[l] <-
            sample(active.terms$sign[as.logical(ass.kegg[active.genes$gene[l], active.terms$term])], 1)
        }
        
        act.genes.trans$sign <- active.genes[act.genes.trans$gene, "sign"]
        act.genes.prot$sign <- active.genes[act.genes.prot$gene, "sign"]
        
        ## Sample background 
        bg <- data.frame(mvrnorm(length(genes.tp), mu = c(0,0), Sigma = sigma.bg)) # from MASS package
        colnames(bg) <- c("bg_T","bg_P")
        rownames(bg) <- genes.tp
        
        # setup data list: init with bg
        sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]] <-
          list(mRNA=bg[genes.t, "bg_T"],
               prot=bg[genes.p, "bg_P"],
               active.terms=active.terms,
               act.terms = act.terms,
               active.genes=active.genes,
               act.genes.trans=act.genes.trans,
               act.genes.prot=act.genes.prot,
               fdr=setNames(c(alpha, beta, sd.bg, sd.sig), c("alpha", "beta","sd.bg","sd.sig")))
        names(sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]][["mRNA"]]) <- genes.t
        names(sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]][["prot"]]) <- genes.p
        
        ## Sample signal 
        #### 3 categories: active in trans, active in prot, active in both
        delta.t <- 1
        delta.p <- 1
        tries <- 0
        while(delta.t>0.03 || delta.p>0.03){
          tries <- tries + 1
          if(tries>=200) break
          #start tries
          sigtp <- intersect(act.genes.trans$gene, act.genes.prot$gene)
          sigt <- act.genes.trans$gene[!(act.genes.trans$gene %in% sigtp)]
          sigp <- act.genes.prot$gene[!(act.genes.prot$gene %in% sigtp)]
          
          signal.tp <- data.frame()
          if(length(sigtp)!=0){
            signal.tp <- data.frame(matrix(mvrnorm(length(sigtp), mu = c(mean.sig, mean.sig), Sigma = sigma.tp), ncol = 2))
            colnames(signal.tp) <- c("sig_T","sig_P")
            rownames(signal.tp) <- sigtp
          }
          
          signal.t <- data.frame(mvrnorm(length(sigt), mu = c(mean.sig, 0), Sigma = sigma.t)) 
          colnames(signal.t) <- c("sig_T","sig_P")
          rownames(signal.t) <- sigt
          
          signal.p <- data.frame(mvrnorm(length(sigp), mu = c(0, mean.sig), Sigma = sigma.p)) 
          colnames(signal.p) <- c("sig_T","sig_P")
          rownames(signal.p) <- sigp
          
          FNs.t <- abs(c(signal.t[sigt, "sig_T"], signal.tp[sigtp, "sig_T"])) < cutoff.a
          delta.t.new <- abs(sum(as.numeric(FNs.t))/length(FNs.t)-beta)
          FNs.p <- abs(c(signal.p[sigp, "sig_P"], signal.tp[sigtp, "sig_P"]))<cutoff.a
          delta.p.new <- abs(sum(as.numeric(FNs.p))/length(FNs.p)-beta)
          
          if((delta.t.new + delta.p.new) < (delta.t + delta.p)){
            # overwrite bg with signal
            if(length(sigtp)!=0){
              sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]][["mRNA"]][intersect(genes.t, sigtp)] <-
                active.genes[intersect(genes.t, sigtp), "sign"]*signal.tp[intersect(genes.t, sigtp), "sig_T"]
              sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]][["prot"]][intersect(genes.p, sigtp)] <-
                active.genes[intersect(genes.p, sigtp), "sign"]*signal.tp[intersect(genes.p, sigtp), "sig_P"]
            }
            
            sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]][["mRNA"]][intersect(genes.t, sigt)] <-
              active.genes[intersect(genes.t, sigt), "sign"]*signal.t[intersect(genes.t, sigt), "sig_T"]
            sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]][["prot"]][intersect(genes.p, sigt)] <-
              active.genes[intersect(genes.p, sigt), "sign"]*signal.t[intersect(genes.p, sigt), "sig_P"]
            
            sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]][["mRNA"]][intersect(genes.t, sigp)] <-
              active.genes[intersect(genes.t, sigp), "sign"]*signal.p[intersect(genes.t, sigp), "sig_T"]
            sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]][["prot"]][intersect(genes.p, sigp)] <-
              active.genes[intersect(genes.p, sigp), "sign"]*signal.p[intersect(genes.p, sigp), "sig_P"]
            
            delta.t <- delta.t.new
            delta.p <- delta.p.new
          }
        }
        
        if(!with_sign){
          sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]]$mRNA[act.genes.trans$gene[act.genes.trans$sign==-1]] <- sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]]$mRNA[act.genes.trans$gene[act.genes.trans$sign==-1]]*-1
          sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]]$prot[act.genes.prot$gene[act.genes.prot$sign==-1]] <- sim.data[[paste0(k)]][[paste0(alpha, "_", beta)]]$prot[act.genes.prot$gene[act.genes.prot$sign==-1]]*-1
        }
      }
    }
  }
  return(sim.data)
}
#' Case Control Simulation
simulate_casecontrols <- function(zscoreSimdata, N=100, case=0.3){
  ### in result matrix: last 30 IDs are cases (for N*case cases). 
  print("Initialising...")
  sim.data <- zscoreSimdata
  ncase <- round(N*case)
  nctrl <- N-ncase
  sim_casectrl <- list()
  print("Starting simulation...")
  print(paste("Estimated time:", length(sim.data)*length(sim.data[[1]])*0.283,"seconds."))
  for(k in (1:length(sim.data))){
    # k=1
    print(paste("Run",k,"of",length(sim.data)))
    level <- sim.data[[k]]
    for(er in (1:length(level))){ # error rate
      # er=1
      cur_sim <- level[[er]]
      a <- 0.1
      b <- 0.1
      sd.case <- 1
      sd.ctrl <- 1
      
      ## Simulate mRNA
      # pro gen in mRNA: wenn akt, mean(ctrl) != mean (case), wenn inaktiv, mean(case) == mean(ctrl)
      # init mit m(case) == m( ctrl)
      .preMRNA <- lapply(cur_sim$mRNA, function(x) rnorm(N, mean = x, sd = sd.ctrl))
      .mRNAm <- matrix( unlist(.preMRNA), length(.preMRNA), byrow=TRUE ) 
      rownames(.mRNAm) <- names(cur_sim$mRNA)
      colnames(.mRNAm) <- c(paste0("ctrl",1:nctrl),paste0("case",1:ncase))
      .preActM <- lapply(cur_sim$mRNA[cur_sim$act.genes.trans$gene], function(x) rnorm(ncase, mean = x+qnorm(1-a, sd=sd.case)+qnorm(1-b, sd=sd.ctrl), sd = sd.ctrl))
      .mRNAm[cur_sim$act.genes.trans$gene,(nctrl+1):N] <- matrix(unlist(.preActM), nrow(cur_sim$act.genes.trans), byrow = T)
      
      ## Simulate prot
      # pro gen in mRNA: wenn akt, mean(ctrl) != mean (case), wenn inaktiv, mean(case) == mean(ctrl)
      # init mit m(case) == m( ctrl)
      .prePROT <- lapply(cur_sim$prot, function(x) rnorm(N, mean = x, sd = sd.ctrl))
      .protm <- matrix( unlist(.prePROT), length(.prePROT), byrow=TRUE )
      rownames(.protm) <- names(cur_sim$prot)
      colnames(.protm) <- c(paste0("ctrl",1:nctrl),paste0("case",1:ncase))
      .preActP <- lapply(cur_sim$prot[cur_sim$act.genes.prot$gene], function(x) rnorm(ncase, mean = x+qnorm(1-a, sd=sd.case)+qnorm(1-b, sd=sd.ctrl), sd = sd.ctrl))
      .protm[cur_sim$act.genes.prot$gene,(nctrl+1):N] <- matrix(unlist(.preActP), nrow(cur_sim$act.genes.prot), byrow = T)
      
      sim_casectrl[[paste0(k)]][[names(level[er])]]$'mRNA' <- .mRNAm
      sim_casectrl[[paste0(k)]][[names(level[er])]]$'prot' <- .protm
      
      .pheno <- cbind(setNames(c(rep(0,nctrl),rep(1,ncase)),c(paste0("ctrl",1:nctrl),paste0("case",1:ncase))),rep(30,100))
      colnames(.pheno) <- c("Case","Age")
      sim_casectrl[[paste0(k)]][[names(level[er])]]$'pheno' <- .pheno
    }
  }
  return(sim_casectrl)
}
#' Plot bivariate zscore distribution of sim.data object
plot_sim_zscores <- function(simdata.z, plotpath){
  ### Plot Trans/Prot Bivariate Distribution
  alpha = 0.1
  beta = 0.1
  for(k in (1:length(simdata.z))){
    # k=1
    level <- simdata.z[[k]]
    for(er in (1:length(level))){ # error rate
      # er=1
      .current_sim <- level[[er]]
      .genes_tp <- intersect(names(.current_sim$mRNA),names(.current_sim$prot))
      .all_tp <- cbind(.current_sim$mRNA[.genes_tp], .current_sim$prot[.genes_tp])
      row.names(.all_tp) <- .genes_tp
      .background_tp <- .all_tp[! row.names(.all_tp) %in% .current_sim$active.genes$gene,]
      .active_tp <- matrix(.all_tp[row.names(.all_tp) %in% intersect(.current_sim$act.genes.trans$gene,.current_sim$act.genes.prot$gene),], ncol = 2)
      .active_t <- .all_tp[row.names(.all_tp) %in% .current_sim$act.genes.trans$gene,]
      .active_p <- .all_tp[row.names(.all_tp) %in% .current_sim$act.genes.prot$gene,]
      colnames(.background_tp) <- c("Trans","Prot")
      colnames(.active_tp) <- c("Trans","Prot")
      colnames(.active_t) <- c("Trans","Prot")
      colnames(.active_p) <- c("Trans","Prot")
      cutoff.plot <- qnorm(1-(.current_sim$fdr[1])/2, sd=.current_sim$fdr[3])
      range.all <- range(c(.background_tp,.active_p,.active_t))
      
      # density
      dir.create(file.path(plotpath), showWarnings = FALSE)
      
      ##
      # ggplot
      amp1 <- data.frame(.background_tp, class="Background")
      if(nrow(.active_tp)!=0) amp2 <- data.frame(.active_tp, class="Active in both")
      amp3 <- data.frame(.active_t, class="Active in Trans")
      amp4 <- data.frame(.active_p, class="Active in Prot")
      amp <- rbind(amp1, amp3, amp4)
      if(nrow(.active_tp)!=0) amp <- rbind(amp, amp2)
      
      .plotle <- ggplot(amp, aes(x=Trans, y=Prot, colour=class)) +
        geom_point() + 
        geom_vline(xintercept = c(-cutoff.plot,cutoff.plot), color="darkgrey") + 
        geom_hline(yintercept = c(-cutoff.plot,cutoff.plot), color="darkgrey") + 
        scale_color_manual(values=c('black','red', 'blue', 'violet'))
      ggsave(paste0(plotpath,.current_sim$fdr[1],"_",.current_sim$fdr[2],"bivariateDens",k,".pdf"), .plotle, width = 16, height = 10, units = "cm")
      ##
      # pdf(paste0(plotpath,.current_sim$fdr[1],"_",.current_sim$fdr[2],"bivariateDens",k,".pdf"))
      # plot(.background_tp, col="black", xlim = range.all, ylim = range.all, main = paste0("alpha = ",.current_sim$fdr[1]," beta = ",.current_sim$fdr[2],"run ",k))
      # abline(v=c(cutoff.plot,-cutoff.plot), h=c(cutoff.plot,-cutoff.plot), col="black")
      # points(.active_t, col="red")
      # points(.active_p, col="blue")
      # points(.active_tp, col="violet")
      #  dev.off() 
    }
  }
}
#' Evaluate: Exact matches
eval_sim <- function(resdata){
  
  fpr <- c(1e-4,c((1:5)/10))
  fnr <- fpr
  
  match = .create_match(fpr)
  
  
  
  
  for(k in 1:length(resdata)){
    # k=1
    
    print(paste("Run",k,"of",length(resdata)))
    for(mrow in 1:nrow(match)){
      # mrow=1
      
      alpha <- match[mrow, "alpha"]
      beta <- match[mrow, "beta"]
      
      ############
      
      ### read Mona results ######
      #test <- readRDS("correlated/mona1trans/1k0.1_0.1.RDS")
      samples <- resdata[[paste0(k)]][[paste0(alpha, "_", beta)]]
      act2D.6 <- sort(samples[order(samples$posterior.mona2d, decreasing = T)[1:6],]$pathway)
      act2D   <- sort(samples[samples$posterior.mona2d>0.6,]$pathway)
      act1Dt  <- sort(samples[samples$posterior.monaMrna>0.6,]$pathway)
      act1Dp  <- sort(samples[samples$posterior.monaProt>0.6,]$pathway)
      act1DOR <- sort(as.character(unique(c(act1Dt, act1Dp))))
      
      actGSEAt <- sort(samples[samples$padj.gseaMrna<0.05,]$pathway)
      actGSEAp <- sort(samples[samples$padj.gseaProt<0.05,]$pathway)
      
      
      act.terms       <- sort(samples[samples$active==1,]$pathway)
      act.terms.trans <- sort(samples[samples$`active.mRNA`==1,]$pathway)
      act.terms.prot  <- sort(samples[samples$`active.prot`==1,]$pathway)
      
      ##############
      
      #### Comparing ########
      if(identical(act2D.6, act.terms)){
        match[match$alpha==alpha & match$beta==beta, "multi.top6"] <-
          match[match$alpha==alpha & match$beta==beta, "multi.top6"] + 1
      }
      if(identical(act2D, act.terms)){
        match[match$alpha==alpha & match$beta==beta, "multi"] <-
          match[match$alpha==alpha & match$beta==beta, "multi"] + 1
      }
      if(identical(act2D, act.terms.trans)){
        match[match$alpha==alpha & match$beta==beta, "multi.trans"] <-
          match[match$alpha==alpha & match$beta==beta, "multi.trans"] + 1
      }
      if(identical(act2D, act.terms.prot)){
        match[match$alpha==alpha & match$beta==beta, "multi.prot"] <-
          match[match$alpha==alpha & match$beta==beta, "multi.prot"] + 1
      }
      if(identical(act1DOR, act.terms)){
        match[match$alpha==alpha & match$beta==beta, "OR"] <-
          match[match$alpha==alpha & match$beta==beta, "OR"] + 1
      }
      if(identical(act1Dt, act.terms)){
        match[match$alpha==alpha & match$beta==beta, "trans"] <-
          match[match$alpha==alpha & match$beta==beta, "trans"] + 1
      }
      if(identical(act1Dt, act.terms.trans)){
        match[match$alpha==alpha & match$beta==beta, "trans.only"] <-
          match[match$alpha==alpha & match$beta==beta, "trans.only"] + 1
      }
      if(identical(act1Dp, act.terms)){
        match[match$alpha==alpha & match$beta==beta, "prot"] <-
          match[match$alpha==alpha & match$beta==beta, "prot"] + 1
      }
      if(identical(act1Dp, act.terms.prot)){
        match[match$alpha==alpha & match$beta==beta, "prot.only"] <-
          match[match$alpha==alpha & match$beta==beta, "prot.only"] + 1
      }
      if(identical(actGSEAt, act.terms.trans)){
        match[match$alpha==alpha & match$beta==beta, "gsea.trans"] <-
          match[match$alpha==alpha & match$beta==beta, "gsea.trans"] + 1
      }
      if(identical(actGSEAp, act.terms.prot)){
        match[match$alpha==alpha & match$beta==beta, "gsea.prot"] <-
          match[match$alpha==alpha & match$beta==beta, "gsea.prot"] + 1
      }
      ###############
    }
    
  }
  return(match)
}
eval_sim2 <- function(resdata){
  
  fpr <- c(1e-4,c((1:5)/10))
  fnr <- fpr
  
  match = .create_match(fpr)
  
  
  
  
  for(k in 1:length(resdata)){
    # k=1
    
    print(paste("Run",k,"of",length(resdata)))
    for(mrow in 1:nrow(match)){
      # mrow=1
      
      alpha <- match[mrow, "alpha"]
      beta <- match[mrow, "beta"]
      
      ############
      
      ### read Mona results ######
      #test <- readRDS("correlated/mona1trans/1k0.1_0.1.RDS")
      samples <- resdata[[paste0(k)]][[paste0(alpha, "_", beta)]]
      act2D.6 <- sort(samples[order(samples$posterior.mona2d, decreasing = T)[1:6],]$pathway)
      act2D   <- sort(samples[samples$posterior.mona2d>0.6,]$pathway)
      # act1Dt  <- sort(samples[samples$posterior.monaMrna>0.6,]$pathway)
      # act1Dp  <- sort(samples[samples$posterior.monaProt>0.6,]$pathway)
      # act1DOR <- sort(as.character(unique(c(act1Dt, act1Dp))))
      
      # actGSEAt <- sort(samples[samples$padj.gseaMrna<0.05,]$pathway)
      # actGSEAp <- sort(samples[samples$padj.gseaProt<0.05,]$pathway)
      
      
      act.terms       <- sort(samples[samples$active==1,]$pathway)
      act.terms.trans <- sort(samples[samples$`active.mRNA`==1,]$pathway)
      act.terms.prot  <- sort(samples[samples$`active.prot`==1,]$pathway)
      
      ##############
      
      #### Comparing ########
      if(identical(act2D.6, act.terms)){
        match[match$alpha==alpha & match$beta==beta, "multi.top6"] <-
          match[match$alpha==alpha & match$beta==beta, "multi.top6"] + 1
      }
      if(identical(act2D, act.terms)){
        match[match$alpha==alpha & match$beta==beta, "multi"] <-
          match[match$alpha==alpha & match$beta==beta, "multi"] + 1
      }
      if(identical(act2D, act.terms.trans)){
        match[match$alpha==alpha & match$beta==beta, "multi.trans"] <-
          match[match$alpha==alpha & match$beta==beta, "multi.trans"] + 1
      }
      if(identical(act2D, act.terms.prot)){
        match[match$alpha==alpha & match$beta==beta, "multi.prot"] <-
          match[match$alpha==alpha & match$beta==beta, "multi.prot"] + 1
      }
      # if(identical(act1DOR, act.terms)){
      #   match[match$alpha==alpha & match$beta==beta, "OR"] <-
      #     match[match$alpha==alpha & match$beta==beta, "OR"] + 1
      # }
      # if(identical(act1Dt, act.terms)){
      #   match[match$alpha==alpha & match$beta==beta, "trans"] <-
      #     match[match$alpha==alpha & match$beta==beta, "trans"] + 1
      # }
      # if(identical(act1Dt, act.terms.trans)){
      #   match[match$alpha==alpha & match$beta==beta, "trans.only"] <-
      #     match[match$alpha==alpha & match$beta==beta, "trans.only"] + 1
      # }
      # if(identical(act1Dp, act.terms)){
      #   match[match$alpha==alpha & match$beta==beta, "prot"] <-
      #     match[match$alpha==alpha & match$beta==beta, "prot"] + 1
      # }
      # if(identical(act1Dp, act.terms.prot)){
      #   match[match$alpha==alpha & match$beta==beta, "prot.only"] <-
      #     match[match$alpha==alpha & match$beta==beta, "prot.only"] + 1
      # }
      # if(identical(actGSEAt, act.terms.trans)){
      #   match[match$alpha==alpha & match$beta==beta, "gsea.trans"] <-
      #     match[match$alpha==alpha & match$beta==beta, "gsea.trans"] + 1
      # }
      # if(identical(actGSEAp, act.terms.prot)){
      #   match[match$alpha==alpha & match$beta==beta, "gsea.prot"] <-
      #     match[match$alpha==alpha & match$beta==beta, "gsea.prot"] + 1
      # }
      ###############
    }
    
  }
  return(match)
}
eval_sim3 <- function(resdata){  
  match = NULL
  for(k in 1:length(resdata)){
    # k=1
    
    print(paste("Run",k,"of",length(resdata)))
    for(mrow in 1:length(resdata[[k]])){
      # mrow=1
      
      # alpha <- match[mrow, "alpha"]
      # beta <- match[mrow, "beta"]
      
      ############
      
      ### read Mona results ######
      #test <- readRDS("correlated/mona1trans/1k0.1_0.1.RDS")
      samples <- resdata[[paste0(k)]][[mrow]]
      # act2D.6 <- sort(samples[order(samples$posterior.mona2d, decreasing = T)[1:6],]$pathway)
      act2D   <- sort(samples[samples$posterior.mona2d>0.6,]$pathway)
      # act1Dt  <- sort(samples[samples$posterior.monaMrna>0.6,]$pathway)
      # act1Dp  <- sort(samples[samples$posterior.monaProt>0.6,]$pathway)
      # act1DOR <- sort(as.character(unique(c(act1Dt, act1Dp))))
      
      # actGSEAt <- sort(samples[samples$padj.gseaMrna<0.05,]$pathway)
      # actGSEAp <- sort(samples[samples$padj.gseaProt<0.05,]$pathway)
      
      
      act.terms       <- sort(samples[samples$active==1,]$pathway)
      act.terms.trans <- sort(samples[samples$`active.mRNA`==1,]$pathway)
      act.terms.prot  <- sort(samples[samples$`active.prot`==1,]$pathway)
      
      ##############
      if(!identical(act2D, act.terms)){
        match <- rbind(match, c(k,0,mrow,sum(as.numeric(samples$posterior.mona2d>0.6))))
      }
      # match <- rbind(match, c(k,1,mrow,sum(as.numeric(samples$posterior.mona2d>0.6))))
    }
    
    
    
    
  }
  return(match)
}
#' Plot Match Matrix
plot_match <- function(matchmat, N){
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
  colorset <- brewer.pal(n = 9, "RdYlBu")
  
  
  a1 <- ggplot(matchmat, aes(x=alpha, y=beta, fill=trans.only/N)) +
    geom_tile() +
    scale_fill_gradient2("exact \nmatch",
                         na.value = "gray",  low = colorset[1],
                         mid=colorset[5],
                         high = colorset[9],midpoint = 0.5,
                         limit = c(0,1)) +
    geom_tile(color='black') +
    theme_bw() +
    xlim(c(-0.1,0.6)) +
    ylim(c(-0.1,0.6)) +
    theme(aspect.ratio = 1) +
    labs(title="Transcriptomics",
         subtitle="(only active trans)") + 
    xlab("FPR") +
    ylab("FNR")
  a2 <- ggplot(matchmat, aes(x=alpha, y=beta, fill=prot.only/N)) +
    geom_tile() +
    scale_fill_gradient2("exact \nmatch",                       # Legende
                         na.value = "gray",  low = colorset[1],
                         mid=colorset[5],
                         high = colorset[9],midpoint = 0.5,
                         limit = c(0,1)) +
    geom_tile(color='black') +                                  # schwarze Gitterstruktur
    theme_bw() +
    xlim(c(-0.1,0.6)) +                                            # Achsen range
    ylim(c(-0.1,0.6)) +
    theme(aspect.ratio = 1) +
    labs(           title="Proteomics",
                    subtitle="(only active prot)") + 
    xlab("FPR") +
    ylab("FNR")
  
  a3 <- ggplot(matchmat, aes(x=alpha, y=beta, fill=trans/N)) +
    geom_tile() +
    scale_fill_gradient2("exact \nmatch",
                         na.value = "gray",  low = colorset[1],
                         mid=colorset[5],
                         high = colorset[9],midpoint = 0.5,
                         limit = c(0,1)) +
    geom_tile(color='black') +
    theme_bw() +
    xlim(c(-0.1,0.6)) +
    ylim(c(-0.1,0.6)) +
    theme(aspect.ratio = 1) +
    labs(title="Transcriptomics",
         subtitle="(all active gene sets)") + 
    xlab("FPR") +
    ylab("FNR")
  a4 <- ggplot(matchmat, aes(x=alpha, y=beta, fill=prot/N)) +
    geom_tile() +
    scale_fill_gradient2("exact \nmatch",
                         na.value = "gray",  low = colorset[1],
                         mid=colorset[5],
                         high = colorset[9],midpoint = 0.5,
                         limit = c(0,1)) +
    geom_tile(color='black') +
    theme_bw() +
    xlim(c(-0.1,0.6)) +
    ylim(c(-0.1,0.6)) +
    theme(aspect.ratio = 1) +
    labs(title="Proteomics",
         subtitle="(all active gene sets)") + 
    xlab("FPR") +
    ylab("FNR")
  
  a5 <- ggplot(matchmat, aes(x=alpha, y=beta, fill=OR/N)) +
    geom_tile() +
    scale_fill_gradient2("exact \nmatch",
                         na.value = "gray",  low = colorset[1],
                         mid=colorset[5],
                         high = colorset[9],midpoint = 0.5,
                         limit = c(0,1)) +
    geom_tile(color='black') +
    theme_bw() +
    xlim(c(-0.1,0.6)) +
    ylim(c(-0.1,0.6)) +
    theme(aspect.ratio = 1) +
    labs( title="Trans- 'OR' Proteomics",
          subtitle="MONA 1D combined") + 
    xlab("FPR") +
    ylab("FNR")
  a6 <- ggplot(matchmat, aes(x=alpha, y=beta, fill=multi/N)) +
    geom_tile() +
    scale_fill_gradient2("exact \nmatch",
                         na.value = "gray",  low = colorset[1],
                         mid=colorset[5],
                         high = colorset[9],midpoint = 0.5,
                         limit = c(0,1)) +
    geom_tile(color='black') +
    theme_bw() +
    xlim(c(-0.1,0.6)) +
    ylim(c(-0.1,0.6)) +
    theme(aspect.ratio = 1) +
    labs(title="Trans- and Proteomics",
         subtitle="MONA 2D") + 
    xlab("FPR") +
    ylab("FNR")
  
  plot.comp <- ggarrange(a1, a2, a3, a4, a5, a6,
                         nrow = 3, ncol = 2,
                         common.legend = T, legend = "right")
}
#########################
### helper functions ####
.ntupel <- function(fpr, fnr=fpr, t=1){
  counter=0
  for (alpha in fpr){
    for(beta in fnr){
      if(alpha+beta>t) break
      counter=counter+1
    }
  }
  return(counter)
}
.create_match <- function(fpr, fnr=fpr){
  .mat <- data.frame(matrix(NA,0,2))
  for (alpha in fpr){
    for(beta in fnr){
      if(alpha+beta>0.7) break
      .mat <- rbind(.mat, c(alpha,beta))
    }
  }
  colnames(.mat) <- c("alpha", "beta")
  .mat$trans <- 0
  .mat$prot <- 0
  .mat$OR <- 0
  .mat$multi <- 0
  .mat$multi.trans <- 0
  .mat$multi.prot <- 0
  .mat$multi.top6 <- 0
  .mat$trans.only <- 0
  .mat$prot.only <- 0
  .mat$gsea.trans <- 0
  .mat$gsea.prot <-0
  return(.mat)
}
.what_cutoff <- function(simdata_level){
  cutoff <- qnorm(1-(simdata_level$fdr[1])/2, sd=simdata_level$fdr[3])
}
.mcFaddens <- function(mod, nullmod){
  1-logLik(mod)/logLik(nullmod)
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
