# setup ------------------------------------------------------------------------
setwd("~/work/multi_omics_enrich/scripts/")

library(ggplot2)
library(ggpubr)
library(plyr)
library(batchtools)

# source files
source("batchtools_helper.R")

iter = 1:100
outdir <- "../results/current/thesis"

pws <- c("shared", "independent")
coverage <- c("real", 0.3, 1)
rho <- c(0.2, 0.3, 0.8)

batch.table <- expand.grid(pws,
                           coverage,
                           rho,
                           stringsAsFactors = F)
batch.table$label <- paste0("cor",
                            "_rho", batch.table$Var3,
                            "_cov", batch.table$Var2,
                            "_", batch.table$Var1)
colnames(batch.table) <- c("effects", "coverage", "correlation", "label")
batch.table$scenario <- batch.table$label

collect_files <- function(kkk){
  
  iter = 1:100
  outdir <- "../results/current/thesis"
  
  pws <- c("shared", "independent")
  coverage <- c("real", 0.3, 1)
  rho <- c(0.2, 0.3, 0.8)
  
  batch.table <- expand.grid(pws,
                             coverage,
                             rho,
                             stringsAsFactors = F)
  batch.table$label <- paste0("cor",
                              "_rho", batch.table$Var3,
                              "_cov", batch.table$Var2,
                              "_", batch.table$Var1)
  colnames(batch.table) <- c("effects", "coverage", "correlation", "label")
  batch.table$scenario <- batch.table$label
  
  for (label in batch.table$label[kkk]){
    # label <- batch.table$label[1]
    
    if(!file.exists(paste0(outdir, "/", label, "/eval/",
                           "evaluation_full.RDS"))){
      print(paste0("Collecting information for file ",
                   label, "/eval/", "evaluation_full.RDS"))
      res <- readRDS(paste0(outdir, "/", label, "/eval/",
                            "evaluation_", 1, ".RDS"))
      for (k in iter[-1]){
        res2 <- readRDS(paste0(outdir, "/", label, "/eval/",
                               "evaluation_", k, ".RDS"))
        res <- rbind(res, res2)
        
        if(k %in% c(10, 20, 30, 40, 50, 60, 70, 80, 90)){
          print(paste0("Iteration ", k, " / 100 done."))
        }
      }
      
      print(paste0("Saving file ",
                   label, "/eval/", "evaluation_full.RDS"))
      saveRDS(res,
              file=paste0(outdir, "/", label, "/eval/",
                          "evaluation_full.RDS"))
      
      print("Aggregate...")
      res.agg <- aggregate(res[, !(colnames(res) %in% c("eval_act_gt", "eval_act_sim", "roc_obj"))],
                           by=list(res$scenario, res$type, res$eval, res$alpha, res$beta),
                           FUN=mean, na.rm=T)
      res.agg[, c("scenario", "type", "eval", "alpha", "beta")] <-
        res.agg[, c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5")]
      res.agg[, c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5", "N")] <- NULL
      
      print(paste0("Saving aggregated results ",
                   label, "/eval/", "evaluation_results.RDS"))
      saveRDS(res.agg,
              file=paste0(outdir, "/", label, "/evaluation_results.RDS"))
      
    }else{
      if(!file.exists(paste0(outdir, "/", label, "/evaluation_results.RDS"))){
        
        print(paste0("Recovering file ",
                     label, "/eval/", "evaluation_full.RDS"))
        res <- readRDS(file=paste0(outdir, "/", label, "/eval/",
                                   "evaluation_full.RDS"))
        
        print("Aggregate...")
        res.agg <- aggregate(res[, !(colnames(res) %in% c("eval_act_gt", "eval_act_sim", "roc_obj"))],
                             by=list(res$scenario, res$type, res$eval, res$alpha, res$beta),
                             FUN=mean, na.rm=T)
        res.agg[, c("scenario", "type", "eval", "alpha", "beta")] <-
          res.agg[, c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5")]
        res.agg[, c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5", "N")] <- NULL
        
        print(paste0("Saving aggregated results ",
                     label, "/eval/", "evaluation_results.RDS"))
        saveRDS(res.agg,
                file=paste0(outdir, "/", label, "/evaluation_results.RDS"))
      }
    }
    print(paste0("Scenario ", label, " done."))
  }
}

library(ggplot2)
library(ggpubr)
library(plyr)

iter = 1:100
outdir <- "../results/current/thesis"

pws <- c("shared", "independent")
coverage <- c("real", 0.3, 1)
rho <- c(0.2, 0.3, 0.8)

batch.table <- expand.grid(pws,
                           coverage,
                           rho,
                           stringsAsFactors = F)
batch.table$label <- paste0("cor",
                            "_rho", batch.table$Var3,
                            "_cov", batch.table$Var2,
                            "_", batch.table$Var1)
colnames(batch.table) <- c("effects", "coverage", "correlation", "label")
batch.table$scenario <- batch.table$label

names <- paste0("sum", 1:dim(batch.table)[1])

run.batchtools(collect_files, 1:dim(batch.table)[1], more.args=list(),
               names, "tmp_dir_eval_thesis_sum",
               clean.up=F,
               resources=list(partition="cpu_p",
                              exclude="cpusrv07,cpusrv26,cpusrv28,ibis216-010-035,ibis216-010-037,ibis216-010-051,ibis216-010-064",
                              memory="120G",
                              ncpus = 1,
                              measure.memory = TRUE,
                              walltime="08:00:00"))

setwd("~/work/multi_omics_enrich/scripts/")
library(batchtools)
loadRegistry(file.dir="tmp_dir_eval_thesis_sum", writeable = F)
getStatus
getJobTable()

if(length(c(findExpired()$job.id, findErrors()$job.id)>0)){
  loadRegistry(file.dir="tmp_dir_eval_thesis_sum", writeable = T)
  getStatus()
  del_ids <- NULL
  del_ids <- c(findExpired()$job.id, findErrors()$job.id)
  if(!is.null(del_ids)){
    submitJobs(del_ids,
               resources=list(partition="cpu_p",
                              exclude="cpusrv07,cpusrv26,cpusrv28,ibis216-010-035,ibis216-010-037,ibis216-010-051,ibis216-010-064",
                              memory="120G",
                              ncpus = 1,
                              measure.memory = TRUE,
                              walltime="08:00:00"))
  }
  getStatus()
}

q(save="no")


# setup ------------------------------------------------------------------------
setwd("~/work/multi_omics_enrich/scripts/")

library(ggplot2)
library(ggpubr)
library(plyr)
library(batchtools)

# source files
source("batchtools_helper.R")

iter = 1:100
outdir <- "../results/current/thesis"

pws <- c("shared", "independent")
coverage <- c("real", 0.3, 1)
rho <- c(0.2, 0.3, 0.8)

batch.table <- expand.grid(pws,
                           coverage,
                           rho,
                           stringsAsFactors = F)
batch.table$label <- paste0("cor",
                            "_rho", batch.table$Var3,
                            "_cov", batch.table$Var2,
                            "_", batch.table$Var1)
colnames(batch.table) <- c("effects", "coverage", "correlation", "label")
batch.table$scenario <- batch.table$label

res.all <- NULL

for(kkk in 1:dim(batch.table)[1]){
  label <- batch.table$label[kkk]
  res.agg <- readRDS(file=paste0(outdir, "/", label, "/evaluation_results.RDS"))
  res.all <- rbind(res.all, res.agg)
}

res.all <- merge(batch.table[, c("effects", "coverage", "correlation", "scenario")],
                 res.all,
                 by = "scenario",
                 all = T)

print("Save final results evaluation_results_all.RDS.")
saveRDS(res.all,
        file=paste0(outdir, "/evaluation_results_all.RDS"))


#res.all <- readRDS(paste0(outdir, "/evaluation_results_all.RDS"))
