# setup ------------------------------------------------------------------------
setwd("~/work/multi_omics_enrich/scripts/")

library(fgsea)
library(ggplot2)
library(ggpubr)
library(batchtools)

# source files
source("batchtools_helper.R")

print(sessionInfo())

run_mona <- function(kkk){
  # for (k in 1:100){
  # kkk = 1
  
  # setup ------------------------------------------------------------------------
  setwd("~/work/multi_omics_enrich/scripts/")
  
  library(fgsea)
  library(ggplot2)
  library(ggpubr)
  library(batchtools)
  library(plyr)
  
  # source files
  source("batchtools_helper.R")
  
  pw_with_dir <- function(Ass){
    A <- Ass
    rownames(A) <- paste0(rownames(Ass), "_up")
    colnames(A) <- paste0(colnames(Ass), "_up")
    
    NAA <- matrix(0, dim(Ass)[1], dim(Ass)[2])
    rownames(NAA) <- paste0(rownames(Ass), "_up")
    colnames(NAA) <- paste0(colnames(Ass), "_down")
    
    NAB <- matrix(0, dim(Ass)[1], dim(Ass)[2])
    rownames(NAB) <- paste0(rownames(Ass), "_down")
    colnames(NAB) <- paste0(colnames(Ass), "_up")
    
    B <- Ass
    rownames(B) <- paste0(rownames(Ass), "_down")
    colnames(B) <- paste0(colnames(Ass), "_down")
    
    Ass_pm <- rbind(cbind(A, NAA), cbind(NAB, B))
    
    return(Ass_pm)
  }
  
  #p.mono <- "cd ~/ext_tools/mono/bin \n ./mono"
  p.mono <- "cd ~/miniconda2/envs/enrich/bin \n ./mono"
  #p.mono <- "./home/icb/ines.assum/ext_tools/mono/bin/mono"
  p.mona <- "~/work/MONA/software/MonaConsoleApp_Andi_mod/MonaConsoleApp/bin/Debug/MonaConsoleApp.exe"
  gmt.file <- "../data/current/gmt_files/c2.cp.kegg.v7.2.symbols.gmt"
  
  ## define simulation settings ---------------------------------------------------
  outdir = "/home/icb/ines.assum/projects/multi_omics_enrich/results/current/thesis"
  print(paste0("Result directory is ", outdir))
  
  pws <- c("shared", "independent")
  coverage <- c("real", 0.3, 1)
  rho <- c(0.2, 0.3, 0.8)
  iter <- 1:100
  
  batch.table <- expand.grid(pws,
                             coverage,
                             rho,
                             iter,
                             stringsAsFactors = F)
  batch.table$names <- paste0("cor",
                              "_rho", batch.table$Var3,
                              "_cov", batch.table$Var2,
                              "_", batch.table$Var1)
  batch.table$job <- paste0(batch.table$Var4,
                            "_cor",
                            "_rho", batch.table$Var3,
                            "_cov", batch.table$Var2,
                            "_", batch.table$Var1)
  
  label = batch.table$names[kkk]
  k <- batch.table$Var4[kkk]
  
  print(paste0("Current scenario is ", label,
               " at iteration ", k))
  
  # define error rates
  fpr <- c(1e-6, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6)
  fnr <- c(1e-6, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6)
  
  error.rates <- data.frame(matrix(NA, nrow = 0, ncol = 5,
                                   dimnames = list(NULL,
                                                   c("id", "alpha", "alpha_rank",
                                                     "beta", "beta_rank"))))
  run=1
  for (i in 1:(length(fpr))){
    # i=1
    alpha <- fpr[i]
    j=0
    while(j<length(fnr)+4-i & j<length(fnr)){
      j <- j+1
      beta <- fnr[j]
      error.rates[run, ] <- c(run, alpha, i, beta, j)
      run <- run+1
    }
  }
  
  # pathway annotations
  gmt <- gmtPathways(gmt.file)
  hidden <- unique(unlist(gmt))
  Ass <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(Ass)[2]){
    Ass[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  sim_all <- readRDS(file=paste0(outdir, "/", label, "/sum_stats/",
                                 "sim_data_all_", k, ".RDS"))
  sim.data <- sim_all[[paste0(error.rates[1, "alpha"], "_", error.rates[1, "beta"])]]
  
  hidden1 <- intersect(sim.data$spec1[!is.na(sim.data$spec1$statistic), "id"], hidden)
  hidden2 <- intersect(sim.data$spec2[!is.na(sim.data$spec2$statistic), "id"], hidden)
  
  Ass1 <- Ass[hidden1, colnames(Ass)[which(colSums(Ass[hidden1,])>=10)]]
  Ass2 <- Ass[hidden2, colnames(Ass)[which(colSums(Ass[hidden2,])>=5)]]
  if(batch.table$Var2[kkk]=="1"){
    Ass2 <- Ass[hidden2, colnames(Ass)[which(colSums(Ass[hidden2,])>=10)]]
  }
  
  terms1 <- colnames(Ass1)
  terms2 <- colnames(Ass2)
  
  ## MONA 1D ---------------------------------------------------------------------
  path <- paste0(outdir, "/", label, "/MONA/MONA_1D_res/", k, "/")
  unlink(path,
         recursive = T)
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
  assign1 <- apply(Ass1, 1, function(x) paste(which(x>0)-1, collapse=","))
  p.assign1 <- paste0(path, "assignmentMatrix_spec1.txt")
  write(assign1, file=p.assign1)
  p.terms1 <- paste0(path, "terms_spec1.txt")
  write.table(terms1, file=p.terms1, col.names=F, row.names = F, quote=F)
  
  assign2 <- apply(Ass2, 1, function(x) paste(which(x>0)-1, collapse=","))
  p.assign2 <- paste0(path, "assignmentMatrix_spec2.txt")
  write(assign2, file=p.assign2)
  p.terms2 <- paste0(path, "terms_spec2.txt")
  write.table(terms2, file=p.terms2, col.names=F, row.names = F, quote=F)

  if(!file.exists(paste0(outdir, "/", label,
                         "/MONA/MONA_1D_res/MONA_1D_res_", k, ".RDS")) |
     file.size(paste0(outdir, "/", label,
                      "/MONA/MONA_1D_res/MONA_1D_res_", k, ".RDS"))<5){
    
    MONA_1D_res <- data.frame()
    for (l in 1:dim(error.rates)[1]){
      #for (l in 1:10){
      # l=1
      alpha <- error.rates[l, "alpha"]
      beta <- error.rates[l, "beta"]
      sim.data <- sim_all[[paste0(alpha, "_", beta)]]
      
      print(paste0("MONA 1D, iteration ", k,
                   ", alpha = ", alpha, ", beta = ", beta,
                   " and active pathways:"))
      print.data.frame(sim.data$meta.data$active.terms, row.names = F)
      
      # MONA 1D species 1
      MONA_temp1 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = names(gmt),
                               "method" = "MONA_1D",
                               stringsAsFactors = F)
      
      spec1 <- as.numeric(sim.data$spec1[hidden1, "significant"])
      p.spec1 <- paste0(path, "spec1.txt")
      write.table(spec1, file=p.spec1,
                  col.names=F, row.names = F, quote=F)
      
      p.out1 <- paste0(path, "output_spec1_", k, "_", alpha, "_", beta, ".txt")
      
      tries <- 0
      while (!file.exists(p.out1) & tries < 11){
        tries <- tries +1
        sys1D1 <- system(paste(p.mono, p.mona, "0",
                              p.assign1, p.spec1, p.terms1, p.out1, "1",
                              sep = " "), intern = T)
      }
      MONAres1 <- read.table(file=p.out1, sep="\t", h=F)
      colnames(MONAres1) <- c("pathway", "posterior_spec1")

      MONA_temp1 <- merge(MONA_temp1, MONAres1, by="pathway", all=T)
      
      # MONA 1D species 2
      MONA_temp2 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = names(gmt),
                               "method" = "MONA_1D",
                               stringsAsFactors = F)
      
      spec2 <- as.numeric(sim.data$spec2[hidden2, "significant"])
      p.spec2 <- paste0(path, "spec2.txt")
      write.table(spec2, file=p.spec2,
                  col.names=F, row.names = F, quote=F)
      
      p.out2 <- paste0(path, "output_spec2_", k, "_", alpha, "_", beta, ".txt")
      
      tries <- 0
      while (!file.exists(p.out2) & tries < 11){
        tries <- tries +1
        sys1D2 <- system(paste(p.mono, p.mona, "0",
                              p.assign2, p.spec2, p.terms2, p.out2, "1",
                              sep = " "), intern = T)
      }
      MONAres2 <- read.table(file=p.out2, sep="\t", h=F)
      colnames(MONAres2) <- c("pathway", "posterior_spec2")
      MONA_temp2 <- merge(MONA_temp2, MONAres2, by="pathway", all=T)
      
      MONA_temp <- merge(MONA_temp1, MONA_temp2,
                         by=c("pathway", "alpha", "beta", "N", "method"),
                         all=T)
      
      print(paste0("MONA 1D, iteration ", k,
                   ", active pathways:"))
      print.data.frame(MONA_temp[which(MONA_temp$posterior_spec1>0.5 |
                                         MONA_temp$posterior_spec2>0.5),
                                 c("pathway", "posterior_spec1", "posterior_spec2")],
                       row.names = F)
      print.data.frame(data.frame(spec1=as.character(data.frame(sys1D1)[c(5, 13, 14, 15, 17), ]),
                                  spec2=as.character(data.frame(sys1D2)[c(5, 13, 14, 15, 17), ])),
                       row.names = F)
      
      MONA_1D_res <- rbind(MONA_1D_res, MONA_temp)
    }
    saveRDS(MONA_1D_res,
            file=paste0(outdir, "/", label, "/MONA/MONA_1D_res/MONA_1D_res_", k, ".RDS"))
  }
  
  unlink(paste0(outdir, "/", label, "/MONA/MONA_1D_res/", k),
         recursive = T)
  
  print(paste0("MONA 1D, iteration ", k, " done for scenario ", label, " (counter = ", kkk, " )."))
  
  
  ## MONA 1Dpm ---------------------------------------------------------------------
  path <- paste0(outdir, "/", label, "/MONA/MONA_1Dpm_res/", k, "/")
  unlink(path,
         recursive = T)
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
  Ass_pm1 <- pw_with_dir(Ass1)
  Ass_pm2 <- pw_with_dir(Ass2)
  
  terms_pm1 <- colnames(Ass_pm1)
  terms_pm2 <- colnames(Ass_pm2)
  
  assign1 <- apply(Ass_pm1, 1, function(x) paste(which(x>0)-1, collapse=","))
  p.assign1 <- paste0(path, "assignmentMatrix_spec1.txt")
  write(assign1, file=p.assign1)
  p.terms1 <- paste0(path, "terms_spec1.txt")
  write.table(terms_pm1, file=p.terms1, col.names=F, row.names = F, quote=F)
  
  assign2 <- apply(Ass_pm2, 1, function(x) paste(which(x>0)-1, collapse=","))
  p.assign2 <- paste0(path, "assignmentMatrix_spec2.txt")
  write(assign2, file=p.assign2)
  p.terms2 <- paste0(path, "terms_spec2.txt")
  write.table(terms_pm2, file=p.terms2, col.names=F, row.names = F, quote=F)
  
  if(!file.exists(paste0(outdir, "/", label,
                         "/MONA/MONA_1Dpm_res/MONA_1Dpm_res_", k, ".RDS")) |
     file.size(paste0(outdir, "/", label,
                      "/MONA/MONA_1Dpm_res/MONA_1Dpm_res_", k, ".RDS"))<5){
    
    MONA_1Dpm_res <- data.frame()
    for (l in 1:dim(error.rates)[1]){
      #for (l in 1:10){
      
      # alpha <- 0.2
      # beta <- 0.1
      # l=1
      alpha <- error.rates[l, "alpha"]
      beta <- error.rates[l, "beta"]
      sim.data <- sim_all[[paste0(alpha, "_", beta)]]
      
      print(paste0("MONA 1Dpm, iteration ", k,
                   ", alpha = ", alpha, ", beta = ", beta,
                   " and active pathways:"))
      print.data.frame(sim.data$meta.data$active.terms, row.names = F)
      
      # MONA 1Dpm species1
      MONA_temp1 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = terms_pm1,
                               "method" = "MONA_1Dpm",
                               stringsAsFactors = F)
      
      spec1 <- c(as.numeric(sim.data$spec1[hidden1, "significant"] &
                              sim.data$spec1[hidden1, "sign"]>0),
                 as.numeric(sim.data$spec1[hidden1, "significant"] &
                              sim.data$spec1[hidden1, "sign"]<0))
      p.spec1 <- paste0(path, "spec1.txt")
      write.table(spec1, file=p.spec1,
                  col.names=F, row.names = F, quote=F)
      
      p.out1 <- paste0(path, "output_spec1_", k, "_", alpha, "_", beta, ".txt")
      
      tries <- 0
      while (!file.exists(p.out1) & tries < 11){
        tries <- tries +1
        sys1Dpm1 <- system(paste(p.mono, p.mona, "0",
                                 p.assign1, p.spec1, p.terms1, p.out1, "1",
                                 sep = " "), intern = T)
      }
      MONAres1 <- read.table(file=p.out1, sep="\t", h=F)
      colnames(MONAres1) <- c("pathway", "posterior_spec1")

      MONA_temp1 <- merge(MONA_temp1, MONAres1, by="pathway", all=T)
      
      # MONA 1Dpm species 2
      MONA_temp2 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = terms_pm2,
                               "method" = "MONA_1Dpm",
                               stringsAsFactors = F)
      
      spec2 <- c(as.numeric(sim.data$spec2[hidden2, "significant"] &
                              sim.data$spec2[hidden2, "sign"]>0),
                 as.numeric(sim.data$spec2[hidden2, "significant"] &
                              sim.data$spec2[hidden2, "sign"]<0))
      p.spec2 <- paste0(path, "spec2.txt")
      write.table(spec2, file=p.spec2,
                  col.names=F, row.names = F, quote=F)
      
      p.out2 <- paste0(path, "output_spec2_", k, "_", alpha, "_", beta, ".txt")
      
      tries <- 0
      while (!file.exists(p.out2) & tries < 11){
        tries <- tries +1
        sys1Dpm2 <- system(paste(p.mono, p.mona, "0",
                                 p.assign2, p.spec2, p.terms2, p.out2, "1",
                                 sep = " "), intern = T)
      }
      MONAres2 <- read.table(file=p.out2, sep="\t", h=F)
      colnames(MONAres2) <- c("pathway", "posterior_spec2")
      MONA_temp2 <- merge(MONA_temp2, MONAres2, by="pathway", all=T)
      
      MONA_temp <- merge(MONA_temp1, MONA_temp2,
                         by=c("pathway", "alpha", "beta", "N", "method"),
                         all=T)
      
      print(paste0("MONA 1Dpm, iteration ", k,
                   ", active pathways:"))
      print.data.frame(MONA_temp[which(MONA_temp$posterior_spec1>0.5 |
                                         MONA_temp$posterior_spec2>0.5),
                                 c("pathway", "posterior_spec1", "posterior_spec2")],
                       row.names = F)
      print.data.frame(data.frame(spec1=as.character(data.frame(sys1Dpm1)[c(5, 13, 14, 15, 17), ]),
                                  spec2=as.character(data.frame(sys1Dpm2)[c(5, 13, 14, 15, 17), ])),
                       row.names = F)
      
      MONA_1Dpm_res <- rbind(MONA_1Dpm_res, MONA_temp)
    }
    saveRDS(MONA_1Dpm_res,
            file=paste0(outdir, "/", label, "/MONA/MONA_1Dpm_res/MONA_1Dpm_res_", k, ".RDS"))
  }
  
  unlink(paste0(outdir, "/", label, "/MONA/MONA_1Dpm_res/", k),
         recursive = T)
  
  print(paste0("MONA 1Dpm, iteration ", k, " done for scenario ", label, " (counter = ", kkk, " )."))

  
  ## MONA 2D ---------------------------------------------------------------------
  path <- paste0(outdir, "/", label, "/MONA/MONA_2D_res/", k, "/")
  unlink(path,
         recursive = T)
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
  assign <- apply(Ass1, 1, function(x) paste(which(x>0)-1, collapse=","))
  p.assign <- paste0(path, "assignmentMatrix_multi.txt")
  write(assign, file=p.assign)
  p.terms <- paste0(path, "terms_multi.txt")
  write.table(terms1, file=p.terms, col.names=F, row.names = F, quote=F)
  
  miss <- c(as.numeric(!(hidden1 %in% hidden2)))
  p.miss <- paste0(path, "missing.txt")
  write.table(miss, file=p.miss, col.names=F, row.names = F, quote=F)
  
  if(!file.exists(paste0(outdir, "/", label,
                         "/MONA/MONA_2D_res/MONA_2D_res_", k, ".RDS")) |
     file.size(paste0(outdir, "/", label,
                      "/MONA/MONA_2D_res/MONA_2D_res_", k, ".RDS"))<5){
    
    MONA_2D_res <- data.frame()
    for (l in 1:dim(error.rates)[1]){
      #for (l in 1:10){
      
      # alpha <- 0.2
      # beta <- 0.1
      # l=1
      alpha <- error.rates[l, "alpha"]
      beta <- error.rates[l, "beta"]
      sim.data <- sim_all[[paste0(alpha, "_", beta)]]
      
      print(paste0("MONA 2D, iteration ", k,
                   ", alpha = ", alpha, ", beta = ", beta,
                   " and active pathways:"))
      print.data.frame(sim.data$meta.data$active.terms, row.names = F)
      
      MONA_temp <- data.frame("alpha" = alpha,
                              "beta" = beta,
                              "N" = k,
                              "pathway" = terms1,
                              "method" = "MONA_2D",
                              stringsAsFactors = F)
      
      spec1 <- as.numeric(sim.data$spec1[hidden1, "significant"])
      spec2 <- as.numeric(sim.data$spec2[hidden1, "significant"])
      
      #test1 <- data.frame(hidden=hidden1, spec1=spec1, spec2=spec2)
      spec2[is.na(spec2)] <- 0
      
      p.spec1 <- paste0(path, "spec1.txt")
      write.table(spec1, file=p.spec1,
                  col.names=F, row.names = F, quote=F)
      p.spec2 <- paste0(path, "spec2.txt")
      write.table(spec2, file=p.spec2,
                  col.names=F, row.names = F, quote=F)
      
      p.out <- paste0(path, "output_multi_", k, "_", alpha, "_", beta, ".txt")
      
      tries <- 0
      while (!file.exists(p.out) & tries < 11){
        tries <- tries +1
        sys2D <- system(paste(p.mono, p.mona, "7",
                              p.assign, p.spec1, p.terms, p.out,
                              "1", p.spec2, p.miss,
                              sep = " "), intern = T)
      }
      MONAres <- read.table(file=p.out, sep="\t", h=F)
      colnames(MONAres) <- c("pathway", "posterior_multi")
      MONA_temp <- merge(MONA_temp, MONAres, by="pathway", all=T)
      
      print(paste0("MONA 2D, iteration ", k,
                   ", active pathways:"))
      print.data.frame(MONA_temp[which(MONA_temp$posterior_multi>0.5),
                                 c("pathway", "posterior_multi")],
                       row.names = F)
      print.data.frame(data.frame(multi=as.character(data.frame(sys2D)[c(7, 19:23, 25), ])),
                       row.names = F)
      
      MONA_2D_res <- rbind(MONA_2D_res, MONA_temp)
    }
    saveRDS(MONA_2D_res,
            file=paste0(outdir, "/", label, "/MONA/MONA_2D_res/MONA_2D_res_", k, ".RDS"))
  }
  
  unlink(paste0(outdir, "/", label, "/MONA/MONA_2D_res/", k),
         recursive = T)
  
  print(paste0("MONA 2D, iteration ", k, " done for scenario ", label, " (counter = ", kkk, " )."))
  
  
  
  ## MONA 2Dpm ---------------------------------------------------------------------
  path <- paste0(outdir, "/", label, "/MONA/MONA_2Dpm_res/", k, "/")
  unlink(path,
         recursive = T)
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
  Ass_pm <- pw_with_dir(Ass1)
  terms_pm <- colnames(Ass_pm)
  
  assign <- apply(Ass_pm, 1, function(x) paste(which(x>0)-1, collapse=","))
  p.assign <- paste0(path, "assignmentMatrix_multi.txt")
  write(assign, file=p.assign)
  p.terms <- paste0(path, "terms_multi.txt")
  write.table(terms_pm, file=p.terms, col.names=F, row.names = F, quote=F)
  
  miss <- c(as.numeric(!(hidden1 %in% hidden2)), as.numeric(!(hidden1 %in% hidden2)))
  p.miss <- paste0(path, "missing.txt")
  write.table(miss, file=p.miss, col.names=F, row.names = F, quote=F)
  
  if(!file.exists(paste0(outdir, "/", label,
                         "/MONA/MONA_2Dpm_res/MONA_2Dpm_res_", k, ".RDS")) |
     file.size(paste0(outdir, "/", label,
                      "/MONA/MONA_2Dpm_res/MONA_2Dpm_res_", k, ".RDS"))<5){
    
    MONA_2Dpm_res <- data.frame()
    for (l in 1:dim(error.rates)[1]){
      #for (l in 1:10){
      
      # alpha <- 0.2
      # beta <- 0.1
      # l=1
      alpha <- error.rates[l, "alpha"]
      beta <- error.rates[l, "beta"]
      sim.data <- sim_all[[paste0(alpha, "_", beta)]]
      
      print(paste0("MONA 2Dpm, iteration ", k,
                   ", alpha = ", alpha, ", beta = ", beta,
                   " and active pathways:"))
      print.data.frame(sim.data$meta.data$active.terms, row.names = F)
      
      MONA_temp <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = terms_pm,
                               "method" = "MONA_2Dpm",
                              stringsAsFactors = F)
      
      spec1 <- c(as.numeric(sim.data$spec1[hidden1, "significant"] &
                              sim.data$spec1[hidden1, "sign"]>0),
                 as.numeric(sim.data$spec1[hidden1, "significant"] &
                              sim.data$spec1[hidden1, "sign"]<0))
      spec2 <- c(as.numeric(sim.data$spec2[hidden1, "significant"] &
                              sim.data$spec2[hidden1, "sign"]>0),
                 as.numeric(sim.data$spec2[hidden1, "significant"] &
                              sim.data$spec2[hidden1, "sign"]<0))
      
      #test1 <- data.frame(hidden=hidden1, spec1=spec1, spec2=spec2)
      spec2[is.na(spec2)] <- 0
      
      p.spec1 <- paste0(path, "spec1.txt")
      write.table(spec1, file=p.spec1,
                  col.names=F, row.names = F, quote=F)
      p.spec2 <- paste0(path, "spec2.txt")
      write.table(spec2, file=p.spec2,
                  col.names=F, row.names = F, quote=F)
      
      p.out <- paste0(path, "output_multi_", k, "_", alpha, "_", beta, ".txt")
      
      tries <- 0
      while (!file.exists(p.out) & tries < 11){
        tries <- tries +1
        sys2Dpm <- system(paste(p.mono, p.mona, "7",
                                p.assign, p.spec1, p.terms, p.out,
                                "1", p.spec2, p.miss,
                                sep = " "), intern = T)
      }
      MONAres <- read.table(file=p.out, sep="\t", h=F)
      colnames(MONAres) <- c("pathway", "posterior_multi")
      MONA_temp <- merge(MONA_temp, MONAres, by="pathway", all=T)
      
      print(paste0("MONA 2Dpm, iteration ", k,
                   ", active pathways:"))
      print.data.frame(MONA_temp[which(MONA_temp$posterior_multi>0.5),
                                 c("pathway", "posterior_multi")],
                       row.names = F)
      print.data.frame(data.frame(multi=as.character(data.frame(sys2Dpm)[c(7, 19:23, 25), ])),
                       row.names = F)
      
      MONA_2Dpm_res <- rbind(MONA_2Dpm_res, MONA_temp)
    }
    saveRDS(MONA_2Dpm_res,
            file=paste0(outdir, "/", label, "/MONA/MONA_2Dpm_res/MONA_2Dpm_res_", k, ".RDS"))
  }
  
  unlink(paste0(outdir, "/", label, "/MONA/MONA_2Dpm_res/", k),
         recursive = T)
  
  print(paste0("MONA 2Dpm, iteration ", k, " done for scenario ", label, " (counter = ", kkk, " )."))
  
  
  print(paste0("MONA complete for iteration ", k, ", scenario ", label, " (counter = ", kkk, " )."))
}

pws <- c("shared", "independent")
coverage <- c("real", 0.3, 1)
rho <- c(0.2, 0.3, 0.8)
iter <- 1:100

batch.table <- expand.grid(pws,
                           coverage,
                           rho,
                           iter,
                           stringsAsFactors = F)
batch.table$names <- paste0("cor",
                            "_rho", batch.table$Var3,
                            "_cov", batch.table$Var2,
                            "_", batch.table$Var1)
batch.table$job <- paste0(batch.table$Var4,
                          "_cor",
                          "_rho", batch.table$Var3,
                          "_cov", batch.table$Var2,
                          "_", batch.table$Var1)
names <- paste0("mona", 1:dim(batch.table)[1])

results <- run.batchtools(run_mona, 1:dim(batch.table)[1], more.args=list(),
                          names, "tmp_dir_mona_thesis",
                          clean.up=F,
                          resources=list(partition="cpu_p",
                                         exclude="ibis216-010-035,cpusrv07,cpusrv11,cpusrv14,cpusrv26,cpusrv28",
                                         memory="1G",
                                         ncpus = 5,
                                         measure.memory = TRUE,
                                         walltime="20:00:00"))

saveRDS(results,
        file = "batchtools_mona_thesis_results.RDS")


loadRegistry(file.dir="tmp_dir_mona_thesis", writeable = T)
del_ids <- NULL
del_ids <- c(del_ids, findExpired()$job.id, findErrors()$job.id)
getStatus()
if(!is.null(del_ids)){
  cbind(batch.table[del_ids, 1:5], No=1:length(del_ids))
  submitJobs(del_ids,
             resources=list(partition="cpu_p",
                            exclude="ibis216-010-035,cpusrv07,cpusrv11,cpusrv14,cpusrv26,cpusrv27,cpusrv28",
                            memory="1G",
                            ncpus = 5,
                            measure.memory = TRUE,
                            walltime="20:00:00"))
}
getStatus()


q(save = "no")

# 
# setwd("~/work/multi_omics_enrich/scripts/")
# 
# library(batchtools)
# 
# loadRegistry(file.dir="tmp_dir_mona_thesis", writeable = T)
# getStatus()
# del_ids <- NULL
# del_ids <- c(findExpired()$job.id, findErrors()$job.id)
# if(!is.null(del_ids)){
#   submitJobs(del_ids,
#              resources=list(partition="cpu_p",
#                             exclude="cpusrv07",
#                             memory="1G",
#                             ncpus = 5,
#                             measure.memory = TRUE,
#                             walltime="20:00:00"))
# }
# getStatus()

# 
# pws <- c("shared", "independent")
# coverage <- c("real", 0.3, 1)
# rho <- c(0.2, 0.3, 0.8)
# iter <- 1:100
# 
# batch.table <- expand.grid(pws,
#                            coverage,
#                            rho,
#                            iter,
#                            stringsAsFactors = F)
# batch.table$names <- paste0("cor",
#                             "_rho", batch.table$Var3,
#                             "_cov", batch.table$Var2,
#                             "_", batch.table$Var1)
# batch.table$job <- paste0(batch.table$Var4,
#                           "_cor",
#                           "_rho", batch.table$Var3,
#                           "_cov", batch.table$Var2,
#                           "_", batch.table$Var1)
# names <- paste0("mona", 1:dim(batch.table)[1])
# 
# library(batchtools)
# setwd("~/work/multi_omics_enrich/scripts/")
# reg <- readRDS(paste0("tmp_dir_mona_thesis", "/registry.rds"))
# findErrors(ids = NULL, reg = reg)
