# setup ------------------------------------------------------------------------
setwd("~/work/multi_omics_enrich/scripts/")

library(fgsea)
library(mgsa)
library(ggplot2)
library(ggpubr)
library(batchtools)

# source files
source("batchtools_helper.R")

print(sessionInfo())

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

#run_mgsa <- function(kkk){
for(kkk in c(35)){  
  # kkk = 2
  
  # setup ------------------------------------------------------------------------
  setwd("~/work/multi_omics_enrich/scripts/")
  
  library(fgsea)
  library(mgsa)
  library(ggplot2)
  library(ggpubr)
  library(batchtools)
  library(plyr)
  
  # source files
  source("batchtools_helper.R")
  
  # Adjacency matrix to list
  matrix_to_list <- function(pws){
    pws.l <- list()
    for (pw in colnames(pws)) {
      pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
    }
    return(pws.l)
  }
  
  # Pathways with direction
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
  
  gmt.file <- "../data/current/gmt_files/c2.cp.kegg.v7.2.symbols.gmt"
  
  ## define simulation settings ---------------------------------------------------
  outdir = "../results/current/thesis"
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
  
  Ass1 <- Ass[hidden1, colnames(Ass)[which(colSums(Ass[hidden1,])>10)]]
  Ass2 <- Ass[hidden2, colnames(Ass)[which(colSums(Ass[hidden2,])>5)]]
  
  terms1 <- colnames(Ass1)
  terms2 <- colnames(Ass2)
  
  pwl1 <- matrix_to_list(Ass1)
  pwl2 <- matrix_to_list(Ass2)
  
  print.data.frame(sim.data$meta.data$active.terms, row.names = F)
  
  ## MGSA P ---------------------------------------------------------------------
  path <- paste0(outdir, "/", label, "/MGSA/MGSA_P_res/")
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }

  if(!file.exists(paste0(outdir, "/", label, "/MGSA/MGSA_P_res/MGSA_P_res_", k, ".RDS")) |
     file.size(paste0(outdir, "/", label, "/MGSA/MGSA_P_res/MGSA_P_res_", k, ".RDS"))<5){
    
    MGSA_P_res <- data.frame()
    MGSA_P_res_list <- list(spec1=list(), spec2=list())
    for (l in 1:dim(error.rates)[1]){
      #for (l in 1:10){
      # l=1
      alpha <- error.rates[l, "alpha"]
      beta <- error.rates[l, "beta"]
      sim.data <- sim_all[[paste0(alpha, "_", beta)]]
      
      print(paste0("MGSA P, iteration ", k,
                   ", alpha = ", alpha, ", beta = ", beta,
                   " and error rates ", l, " of ", dim(error.rates)[1]))
      
      # MGSA P species 1
      MGSA_temp1 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = names(gmt),
                               "method" = "MGSA_P",
                               stringsAsFactors = F)
      
      spec1 <- sim.data$spec1$id[which(sim.data$spec1$significant==1)]
      
      res <- NULL
      tries <- 0
      while(is.null(res) & tries < 10){
        tries <- tries + 1
        try({res <- mgsa(spec1, pwl1)})
      }
      
      MGSA_P_res_list[["spec1"]][[paste0(alpha, "_", beta)]] <- res
      MGSAres1 <- res@setsResults
      MGSAres1$pathway <- rownames(MGSAres1)
      colnames(MGSAres1) <- c("inPopulation_spec1", "inStudySet_spec1", "posterior_spec1", "std.error_spec1", "pathway")
      
      MGSA_temp1 <- merge(MGSA_temp1, MGSAres1, by="pathway", all=T)
      
      # MGSA P species 2
      MGSA_temp2 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = names(gmt),
                               "method" = "MGSA_P",
                               stringsAsFactors = F)
      
      spec2 <- sim.data$spec2$id[which(sim.data$spec2$significant==1)]
      
      res <- NULL
      tries <- 0
      while(is.null(res) & tries < 10){
        tries <- tries + 1
        try({res <- mgsa(spec2, pwl2)})
      }

      MGSA_P_res_list[["spec2"]][[paste0(alpha, "_", beta)]] <- res
      MGSAres2 <- res@setsResults
      MGSAres2$pathway <- rownames(MGSAres2)
      colnames(MGSAres2) <- c("inPopulation_spec2", "inStudySet_spec2", "posterior_spec2", "std.error_spec2", "pathway")
      
      MGSA_temp2 <- merge(MGSA_temp2, MGSAres2, by="pathway", all=T)
      
      MGSA_temp <- merge(MGSA_temp1, MGSA_temp2,
                         by=c("pathway", "alpha", "beta", "N", "method"),
                         all=T)
      
      print(paste0("MGSA P, iteration ", k,
                   ", active pathways:"))
      print.data.frame(MGSA_temp[which(MGSA_temp$posterior_spec1>0.5 |
                                         MGSA_temp$posterior_spec2>0.5),
                                 c("pathway", "posterior_spec1", "posterior_spec2")],
                       row.names = F)
      
      MGSA_P_res <- rbind(MGSA_P_res, MGSA_temp)
    }
    saveRDS(MGSA_P_res,
            file=paste0(outdir, "/", label, "/MGSA/MGSA_P_res/MGSA_P_res_", k, ".RDS"))
    saveRDS(MGSA_P_res_list,
            file=paste0(outdir, "/", label, "/MGSA/MGSA_P_res/MGSA_P_res_list_", k, ".RDS"))
  }
  
  print(paste0("MGSA P, iteration ", k, " done for scenario ", label, " (counter = ", kkk, " )."))
  
  
  ## MGSA T ---------------------------------------------------------------------
  path <- paste0(outdir, "/", label, "/MGSA/MGSA_T_res/")
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
  Ass_pm1 <- pw_with_dir(Ass1)
  Ass_pm2 <- pw_with_dir(Ass2)
  
  terms_pm1 <- colnames(Ass_pm1)
  terms_pm2 <- colnames(Ass_pm2)
  
  pwl_pm1 <- matrix_to_list(Ass_pm1)
  pwl_pm2 <- matrix_to_list(Ass_pm2)
  
  if(!file.exists(paste0(outdir, "/", label, "/MGSA/MGSA_T_res/MGSA_T_res_", k, ".RDS")) |
     file.size(paste0(outdir, "/", label, "/MGSA/MGSA_T_res/MGSA_T_res_", k, ".RDS"))<5){
    
    MGSA_T_res <- data.frame()
    MGSA_T_res_list <- list(spec1=list(), spec2=list())
    for (l in 1:dim(error.rates)[1]){
      #for (l in 1:10){
      
      # alpha <- 0.2
      # beta <- 0.1
      # l=1
      alpha <- error.rates[l, "alpha"]
      beta <- error.rates[l, "beta"]
      sim.data <- sim_all[[paste0(alpha, "_", beta)]]
      
      print(paste0("MGSA T, iteration ", k,
                   ", alpha = ", alpha, ", beta = ", beta,
                   " and error rates ", l, " of ", dim(error.rates)[1]))
      
      # MGSA T species1
      MGSA_temp1 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = terms_pm1,
                               "method" = "MGSA_T",
                               stringsAsFactors = F)
      
      spec1 <- paste0(sim.data$spec1$id[which(sim.data$spec1$significant==1)],
                      ifelse(sim.data$spec1$sign[which(sim.data$spec1$significant==1)]>0,
                             "_up", "_down"))
      
      res <- NULL
      tries <- 0
      while(is.null(res) & tries < 10){
        tries <- tries + 1
        try({res <- mgsa(spec1, pwl_pm1)})
      }
      
      MGSA_T_res_list[["spec1"]][[paste0(alpha, "_", beta)]] <- res
      MGSAres1 <- res@setsResults
      MGSAres1$pathway <- rownames(MGSAres1)
      colnames(MGSAres1) <- c("inPopulation_spec1", "inStudySet_spec1",
                              "posterior_spec1", "std.error_spec1",
                              "pathway")
      
      MGSA_temp1 <- merge(MGSA_temp1, MGSAres1, by="pathway", all=T)
      
      # MGSA T species 2
      MGSA_temp2 <- data.frame("alpha" = alpha,
                               "beta" = beta,
                               "N" = k,
                               "pathway" = terms_pm2,
                               "method" = "MGSA_T",
                               stringsAsFactors = F)
      
      spec2 <- paste0(sim.data$spec2$id[which(sim.data$spec2$significant==1)],
                      ifelse(sim.data$spec2$sign[which(sim.data$spec2$significant==1)]>0,
                             "_up", "_down"))
      
      res <- NULL
      tries <- 0
      while(is.null(res) & tries < 10){
        tries <- tries + 1
        try({res <- mgsa(spec2, pwl_pm2)})
      }
      MGSA_T_res_list[["spec2"]][[paste0(alpha, "_", beta)]] <- res
      MGSAres2 <- res@setsResults
      MGSAres2$pathway <- rownames(MGSAres2)
      colnames(MGSAres2) <- c("inPopulation_spec2", "inStudySet_spec2",
                              "posterior_spec2", "std.error_spec2",
                              "pathway")
      
      MGSA_temp2 <- merge(MGSA_temp2, MGSAres2, by="pathway", all=T)
      
      MGSA_temp <- merge(MGSA_temp1, MGSA_temp2,
                         by=c("pathway", "alpha", "beta", "N", "method"),
                         all=T)
      
      print(paste0("MGSA T, iteration ", k,
                   ", active pathways:"))
      print.data.frame(MGSA_temp[which(MGSA_temp$posterior_spec1>0.5 |
                                         MGSA_temp$posterior_spec2>0.5),
                                 c("pathway", "posterior_spec1", "posterior_spec2")],
                       row.names = F)
      
      MGSA_T_res <- rbind(MGSA_T_res, MGSA_temp)
    }
    saveRDS(MGSA_T_res,
            file=paste0(outdir, "/", label, "/MGSA/MGSA_T_res/MGSA_T_res_", k, ".RDS"))
    saveRDS(MGSA_T_res_list,
            file=paste0(outdir, "/", label, "/MGSA/MGSA_T_res/MGSA_T_res_list_", k, ".RDS"))
  }
  
  print(paste0("MGSA T, iteration ", k, " done for scenario ", label, " (counter = ", kkk, " )."))
  
  
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
names <- paste0("mgsa", 1:dim(batch.table)[1])

results <- run.batchtools(run_mgsa, 1:dim(batch.table)[1], more.args=list(),
                          names, "tmp_dir_mgsa_thesis",
                          clean.up=F,
                          resources=list(partition="cpu_p",
                                         exclude="cpusrv07",
                                         memory="500M",
                                         ncpus = 5,
                                         measure.memory = TRUE,
                                         walltime="01:00:00"))

saveRDS(results,
        file = "batchtools_mgsa_thesis_results.RDS")


loadRegistry(file.dir="tmp_dir_mgsa_thesis", writeable = T)
del_ids <- c(findExpired()$job.id, findError()$job.id)
if(!is.null(del_ids)){
  cbind(batch.table[del_ids, 1:5], No=1:length(del_ids))
  submitJobs(del_ids,
             resources=list(partition="cpu_p",
                            exclude="cpusrv07",
                            memory="500M",
                            ncpus = 5,
                            measure.memory = TRUE,
                            walltime="01:00:00"))
}



q(save = "no")




setwd("~/work/multi_omics_enrich/scripts/")

library(batchtools)

loadRegistry(file.dir="tmp_dir_mgsa_thesis", writeable = T)
getStatus()
del_ids <- NULL
del_ids <- c(findExpired()$job.id, findErrors()$job.id)
if(!is.null(del_ids)){
  submitJobs(del_ids,
             resources=list(partition="cpu_p",
                            exclude="cpusrv07",
                            memory="500M",
                            ncpus = 5,
                            measure.memory = TRUE,
                            walltime="01:00:00"))
}
getStatus()

