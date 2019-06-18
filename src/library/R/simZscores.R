#!/usr/bin/env Rscript
library(fgsea)
library(MASS)
library(pROC)
# library(svMisc)
### OptParse ########
library(optparse)
option_list <- list( 
  make_option(c( "--out1"), type="character", default="mRNA.rds", 
              help="Output file name for mRNA. (Default: hsa"), 
  make_option( "--out2", type="character", default="prot.rds",
              help="Output file name for proteins. (Default: kegg)"),
  make_option(c("--outT"), type="character", default="prot.rds",
              help="Output file name for ground truth. (default= groundtruth.rds)"),
  make_option(c("--gmt"), type="character", default="/data/gmt/c2.cp.kegg.v6.0.symbols.gmt", 
              help="Path to gmt file. (Default: /data/gmt/c2.cp.kegg.v6.0.symbols.gmt"), 
  make_option( "--rho", type="double", default="0.5",
               help="Correalation coeff of species. (Default: 0.5)"),
  make_option(c("--cov"), type="double", default="0.25",
              help="First species coverage. (default= 0.25)"),
  make_option(c( "--min"), type="double", default="10", 
              help="minSize of pathways. (Default: 10"), 
  make_option( "--alpha", type="double", default="kegg",
               help="Error rate alpha. (Default: 0.1)"),
  make_option(c( "--beta"), type="double", default="pathways.gmt",
              help="Error rate beta. (Default: 0.1)"),
  make_option(c("--sign"), type="character", default="yes", 
              help="Include downregulated pathways? (Default: yes)"), 
  make_option( "--sigma", type="double", default="10",
               help="Standard deviation sigma of signal distribution. (Default: 10)")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


out1     <- opt$out1
out2     <- opt$out2
outT     <- opt$outT
gmtpath <- opt$gmt
rho <- opt$rho
coverage <- opt$cov
minsize <- opt$min
alpha <- opt$alpha
beta <- opt$beta
include_direction <- if(opt$sign=="yes")
sdsig <- opt$sigma
#####################
### Functions #######
simulate_single_zscores <- function(gmtpath, rho=0.5, coverage=0.25, minsize=10, alpha, beta, include_direction=T, dynamic_SD=F, sdsig=5){
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
  
  
  genelist <- readRDS("/data/example_data/gene_list_new.RDS")
  genes.t <- intersect(genelist$trans, genes.kegg)
  genes.p <- intersect(genelist$prot25, genes.kegg)
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
  sim.data <-
    list(mRNA=bg[genes.t, "bg_T"],
         prot=bg[genes.p, "bg_P"],
         active.terms=active.terms,
         act.terms = act.terms,
         active.genes=active.genes,
         act.genes.trans=act.genes.trans,
         act.genes.prot=act.genes.prot,
         fdr=setNames(c(alpha, beta, sd.bg, sd.sig), c("alpha", "beta","sd.bg","sd.sig")))
  names(sim.data[["mRNA"]]) <- genes.t
  names(sim.data[["prot"]]) <- genes.p
  
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
        sim.data[["mRNA"]][intersect(genes.t, sigtp)] <-
          active.genes[intersect(genes.t, sigtp), "sign"]*signal.tp[intersect(genes.t, sigtp), "sig_T"]
        sim.data[["prot"]][intersect(genes.p, sigtp)] <-
          active.genes[intersect(genes.p, sigtp), "sign"]*signal.tp[intersect(genes.p, sigtp), "sig_P"]
      }
      
      sim.data[["mRNA"]][intersect(genes.t, sigt)] <-
        active.genes[intersect(genes.t, sigt), "sign"]*signal.t[intersect(genes.t, sigt), "sig_T"]
      sim.data[["prot"]][intersect(genes.p, sigt)] <-
        active.genes[intersect(genes.p, sigt), "sign"]*signal.t[intersect(genes.p, sigt), "sig_P"]
      
      sim.data[["mRNA"]][intersect(genes.t, sigp)] <-
        active.genes[intersect(genes.t, sigp), "sign"]*signal.p[intersect(genes.t, sigp), "sig_T"]
      sim.data[["prot"]][intersect(genes.p, sigp)] <-
        active.genes[intersect(genes.p, sigp), "sign"]*signal.p[intersect(genes.p, sigp), "sig_P"]
      
      delta.t <- delta.t.new
      delta.p <- delta.p.new
    }
  }
  
  if(!with_sign){
    sim.data$mRNA[act.genes.trans$gene[act.genes.trans$sign==-1]] <- sim.data$mRNA[act.genes.trans$gene[act.genes.trans$sign==-1]]*-1
    sim.data$prot[act.genes.prot$gene[act.genes.prot$sign==-1]] <- sim.data$prot[act.genes.prot$gene[act.genes.prot$sign==-1]]*-1
  }
  
  
  return(sim.data)
}
#####################
### Main ############
print(paste("Simulating..."))
result <- simulate_single_zscores(gmtpath, rho, coverage, minsize, alpha, beta, include_direction, dynamic_SD = F,sdsig)
saveRDS(result$mRNA, file = out1)
saveRDS(result$prot, file = out2)
saveRDS(result, file = outT)
print(paste("Done!"))
#####################
# ## Debug
# out1     <- "mRNA.rds"
# out2     <- "prot.rds"
# outT     <- "groundtruth.rds"
# gmtpath <- "/Users/kris.g/maConstruction/scripts/SimEval/c2.cp.kegg.v6.0.symbols.gmt"
# rho <- 0.5
# coverage <- 0.25
# minsize <- 10
# alpha <- 0.1
# beta <- 0.1
# include_direction <- T
# sdsig <- 10
