# setup ------------------------------------------------------------------------
setwd("~/work/multi_omics_enrich/scripts/")

library(fgsea)
library(ggplot2)
library(ggpubr)
library(batchtools)

# source files
source("batchtools_helper.R")

print(sessionInfo())

run_gsea <- function(kkk){
  # for (k in 1:100){
  # kkk = 1
  # kkk = 1488
  
  # setup ------------------------------------------------------------------------
  setwd("~/work/multi_omics_enrich/scripts/")
  
  library(fgsea)
  library(ggplot2)
  library(ggpubr)
  library(batchtools)
  library(plyr)
  
  # source files
  source("batchtools_helper.R")
  
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
  
  
  ## define simulation settings ---------------------------------------------------
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
  
  label = batch.table$names[kkk]
  k <- batch.table$Var4[kkk]
  outdir = "../results/current/thesis"
  
  if (!dir.exists(paste0(outdir, "/", label, "/GSEA/GSEA_P_res/"))) {
    dir.create(paste0(outdir, "/", label, "/GSEA/GSEA_P_res/"), recursive = T)
  }
  
  if(!file.exists(paste0(outdir, "/", label, "/GSEA/GSEA_P_res/GSEA_P_res_", k, ".RDS")) |
     file.size(paste0(outdir, "/", label, "/GSEA/GSEA_P_res/GSEA_P_res_", k, ".RDS"))<5){
    # get pathway file
    
    #test10 <- system.time({
    gmt <- gmtPathways("../data/current/gmt_files/c2.cp.kegg.v7.2.symbols.gmt")
    sim_all <- readRDS(file=paste0(outdir, "/", label, "/sum_stats/",
                                   "sim_data_all_", k, ".RDS"))
    
    GSEA_P_res <- data.frame()
    for (l in 1:dim(error.rates)[1]){
      #for (l in 1:10){
      
      # alpha <- 0.2
      # beta <- 0.1
      # l=1
      alpha <- error.rates[l, "alpha"]
      beta <- error.rates[l, "beta"]
      sim.data <- sim_all[[paste0(alpha, "_", beta)]]
      
      print(paste0("GSEA P, iteration ", k,
                   ", alpha = ", alpha, ", beta = ", beta,
                   " and error rates ", l, " of ", dim(error.rates)[1]))
      
      GSEA_temp1 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = names(gmt),
                               "method" = "GSEA_P")
      
      rank <- sim.data[["spec1"]][complete.cases(sim.data[["spec1"]]$statistic), ]
      rank1 <- abs(rank$statistic)
      names(rank1) <- rank$id
      
      GSEAres1 <- NULL
      GSEAres1 <- fgsea(pathways = gmt,
                        stats = rank1, scoreType = "pos",
                        minSize=10,
                        maxSize=500,
                        eps = 0,
                        nproc = 10)
      
      GSEAres1 <- rename(GSEAres1, c("pval" = "pval_spec1", "padj" = "padj_spec1",
                                     "log2err" = "log2err_spec1",
                                     "ES" = "ES_spec1", "NES" = "NES_spec1",
                                     #"nMoreExtreme" = "nMoreExtreme_spec1",
                                     "size" = "size_spec1",
                                     "leadingEdge" = "leadingEdge_spec1"))
      GSEA_temp1 <- merge(GSEA_temp1, GSEAres1, by="pathway", all=T)
      
      GSEA_temp2 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = names(gmt),
                               "method" = "GSEA_P")
      
      rank <- sim.data[["spec2"]][complete.cases(sim.data[["spec2"]]$statistic), ]
      rank2 <- abs(rank$statistic)
      names(rank2) <- rank$id
      
      GSEAres2 <- NULL
      minsize <- ifelse(batch.table$Var2[kkk]=="1", 10, 5)
      GSEAres2 <- fgsea(pathways = gmt,
                        stats = rank2, scoreType = "pos",
                        minSize=minsize,
                        maxSize=500,
                        eps = 0,
                        nproc = 10)
      
      GSEAres2 <- rename(GSEAres2, c("pval" = "pval_spec2", "padj" = "padj_spec2",
                                     "log2err" = "log2err_spec2",
                                     "ES" = "ES_spec2", "NES" = "NES_spec2",
                                     #"nMoreExtreme" = "nMoreExtreme_spec2",
                                     "size" = "size_spec2",
                                     "leadingEdge" = "leadingEdge_spec2"))
      GSEA_temp2 <- merge(GSEA_temp2, GSEAres2, by="pathway", all=T)
      GSEA_temp <- merge(GSEA_temp1, GSEA_temp2,
                         by=c("pathway", "alpha", "beta", "N", "method"),
                         all=T)
      
      GSEA_P_res <- rbind(GSEA_P_res, GSEA_temp)
    }
    #})
    saveRDS(GSEA_P_res,
            file=paste0(outdir, "/", label, "/GSEA/GSEA_P_res/GSEA_P_res_", k, ".RDS"))
  }

  print(paste0("GSEA P value, iteration ", k, " done for scenario ", label, " (counter = ", kkk, " )."))
  
  if (!dir.exists(paste0(outdir, "/", label, "/GSEA/GSEA_T_res/"))) {
    dir.create(paste0(outdir, "/", label, "/GSEA/GSEA_T_res/"), recursive = T)
  }
  
  if(!file.exists(paste0(outdir, "/", label, "/GSEA/GSEA_T_res/GSEA_T_res_", k, ".RDS")) |
     file.size(paste0(outdir, "/", label, "/GSEA/GSEA_T_res/GSEA_T_res_", k, ".RDS"))<5){
    # get pathway file
    gmt <- gmtPathways("../data/current/gmt_files/c2.cp.kegg.v7.2.symbols.gmt")
    sim_all <- readRDS(file=paste0(outdir, "/", label, "/sum_stats/",
                                   "sim_data_all_", k, ".RDS"))
    
    GSEA_T_res <- data.frame()
    for (l in 1:dim(error.rates)[1]){
      #for (l in 1:10){
      
      # alpha <- 0.2
      # beta <- 0.1
      # l=1
      alpha <- error.rates[l, "alpha"]
      beta <- error.rates[l, "beta"]
      sim.data <- sim_all[[paste0(alpha, "_", beta)]]
      
      print(paste0("GSEA T, iteration ", k,
                   ", alpha = ", alpha, ", beta = ", beta,
                   " and error rates ", l, " of ", dim(error.rates)[1]))
      
      GSEA_temp1 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = names(gmt),
                               "method" = "GSEA_T")
      
      rank <- sim.data[["spec1"]][complete.cases(sim.data[["spec1"]]$statistic), ]
      rank1 <- rank$statistic
      names(rank1) <- rank$id
      
      GSEAres1 <- NULL
      GSEAres1 <- fgsea(pathways = gmt,
                        stats = rank1,
                        minSize = 10,
                        maxSize = 500,
                        eps = 0,
                        nproc = 10)
      
      GSEAres1 <- rename(GSEAres1, c("pval" = "pval_spec1", "padj" = "padj_spec1",
                                     "log2err" = "log2err_spec1",
                                     "ES" = "ES_spec1", "NES" = "NES_spec1",
                                     #"nMoreExtreme" = "nMoreExtreme_spec1",
                                     "size" = "size_spec1",
                                     "leadingEdge" = "leadingEdge_spec1"))
      GSEA_temp1 <- merge(GSEA_temp1, GSEAres1, by="pathway", all=T)
      
      GSEA_temp2 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = names(gmt),
                               "method" = "GSEA_T")
      
      rank <- sim.data[["spec2"]][complete.cases(sim.data[["spec2"]]$statistic), ]
      rank2 <- rank$statistic
      names(rank2) <- rank$id
      
      GSEAres2 <- NULL
      minsize <- ifelse(batch.table$Var2[kkk]=="1", 10, 5)
      GSEAres2 <- fgsea(pathways = gmt,
                        stats = rank2,
                        minSize=minsize,
                        maxSize=500,
                        eps = 0,
                        nproc = 10)
      
      GSEAres2 <- rename(GSEAres2, c("pval" = "pval_spec2", "padj" = "padj_spec2",
                                     "log2err" = "log2err_spec2",
                                     "ES" = "ES_spec2", "NES" = "NES_spec2",
                                     #"nMoreExtreme" = "nMoreExtreme_spec2",
                                     "size" = "size_spec2",
                                     "leadingEdge" = "leadingEdge_spec2"))
      GSEA_temp2 <- merge(GSEA_temp2, GSEAres2, by="pathway", all=T)
      
      GSEA_temp <- merge(GSEA_temp1, GSEA_temp2,
                         by=c("pathway", "alpha", "beta", "N", "method"),
                         all=T)
      
      GSEA_T_res <- rbind(GSEA_T_res, GSEA_temp)
    }
    saveRDS(GSEA_T_res,
            file=paste0(outdir, "/", label, "/GSEA/GSEA_T_res/GSEA_T_res_", k, ".RDS"))
  }
  print(paste0("GSEA T value, iteration ", k, " done for scenario ", label, " (counter = ", kkk, " )."))
  
  print(paste0("GSEA complete for iteration ", k, ", scenario ", label, " (counter = ", kkk, " )."))
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
names <- paste0("gsea", 1:dim(batch.table)[1])

results <- run.batchtools(run_gsea, 1:dim(batch.table)[1], more.args=list(),
                          names, "tmp_dir_gsea_thesis",
                          clean.up=F,
                          resources=list(partition="cpu_p",
                                         exclude = "ibis216-010-035,ibis216-010-064,cpusrv07,cpusrv26,cpusrv27,cpusrv28",
                                         memory="1G",
                                         ncpus = 10,
                                         measure.memory = TRUE,
                                         walltime="04:00:00"))

saveRDS(results,
        file = "batchtools_gsea_thesis_results.RDS")


# Redo expires -----------------------------------------

setwd("~/work/multi_omics_enrich/scripts/")
library(batchtools)
outdir = "../results/current/thesis"
loadRegistry(file.dir="tmp_dir_gsea_thesis", writeable = F)

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

loadRegistry(file.dir="tmp_dir_gsea_thesis", writeable = T)
del_ids <- NULL
del_ids <- unique(c(del_ids, findExpired()$job.id, findErrors()$job.id))
if(!is.null(del_ids)){
  #cbind(batch.table[del_ids, 1:5], No=1:length(del_ids))
  submitJobs(del_ids,
             resources=list(partition="cpu_p",
                            exclude = "ibis216-010-035,ibis216-010-064,cpusrv07,cpusrv26,cpusrv28",
                            memory="1G",
                            ncpus = 10,
                            measure.memory = TRUE,
                            walltime="04:00:00"))
}

q(save = "no")
