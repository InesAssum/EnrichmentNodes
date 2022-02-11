# setup ------------------------------------------------------------------------
setwd("~/work/multi_omics_enrich/scripts/")

library(fgsea)
library(ggplot2)
library(ggpubr)
library(batchtools)

# source files
source("simulate_summary_statistics.R")
source("batchtools_helper.R")

print(sessionInfo())

simulate_data <- function(kkk){
  # setup ------------------------------------------------------------------------
  setwd("~/work/multi_omics_enrich/scripts/")
  
  library(fgsea)
  library(ggplot2)
  library(ggpubr)
  library(batchtools)
  
  # source files
  source("simulate_summary_statistics.R")
  source("batchtools_helper.R")
  
  pws <- c("shared", "independent")
  coverage <- c("real", 0.3, 1)
  rho <- c(0.2, 0.3, 0.8)
  
  batch.table <- expand.grid(pws,
                             coverage,
                             rho,
                             stringsAsFactors = F)
  batch.table$names <- paste0("cor",
                              "_rho", batch.table$Var3,
                              "_cov", batch.table$Var2,
                              "_", batch.table$Var1)
  
  ## define simulation settings ---------------------------------------------------
  
  # get pathway file
  gmt <- gmtPathways("../data/current/gmt_files/c2.cp.kegg.v7.2.symbols.gmt")
  hidden <- unique(unlist(gmt))
  Ass <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(Ass)[2]){
    Ass[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
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
  
  # set parameters
  if(batch.table$Var2[kkk]=="real"){
    cov <- "real"
    cov.path <- paste0("~/work/multi_omics_enrich/data/current/",
                       "coverage_data/coverage_AFHRI_data.RDS")
  }else{
    cov <- as.numeric(batch.table$Var2[kkk])
    cov.path <- NA
  }
  
  if(batch.table$Var1[kkk]=="independent"){
    scen <- "independent"
    nact <- list(spec1=3, spec2=3)
  }else if(batch.table$Var1[kkk]=="shared"){
    scen <- "shared"
    nact <- list(shared=6)
  }
  
  sdbg = 2
  sdsig = "1.5-alpha"
  
  errors1 <- simulate_summary_stats(label=batch.table$names[kkk],
                                    outdir="../results/current/thesis",
                                    seed=NA,
                                    iterations=1:100, error.rates,
                                    gmt,
                                    gmt.path="../data/current/gmt_files/c2.cp.kegg.v7.2.symbols.gmt",
                                    coverage = cov,
                                    coverage.path = cov.path,
                                    scenario=scen,
                                    nactive = nact,
                                    rho=batch.table$Var3[kkk],
                                    filter.miss=10,
                                    sd.bg = sdbg,
                                    sd.signal = sdsig)
}

pws <- c("shared", "independent")
coverage <- c("real", 0.3, 1)
rho <- c(0.2, 0.3, 0.8)

batch.table <- expand.grid(pws,
                           coverage,
                           rho,
                           stringsAsFactors = F)
batch.table$names <- paste0("cor",
                            "_rho", batch.table$Var3,
                            "_cov", batch.table$Var2,
                            "_", batch.table$Var1)
names <- paste0("SimSum", 1:dim(batch.table)[1])

run.batchtools(simulate_data, 1:dim(batch.table)[1], more.args=list(),
               names, "tmp_dir_SimSum",
               clean.up=F,
               resources=list(partition="cpu_p",
                              memory="500M",
                              ncpus = 5,
                              measure.memory = TRUE,
                              walltime="48:00:00"))





