setwd("~/work/multi_omics_enrich/scripts/")

library(fgsea)
library(pROC)
library(ggplot2)
library(ggpubr)
library(batchtools)

# source files
source("batchtools_helper.R")

print(sessionInfo())

run_eval <- function(kkk){
  # for (k in 1:100){
  # kkk = 1
  
  # setup ------------------------------------------------------------------------
  setwd("~/work/multi_omics_enrich/scripts/")
  
  library(fgsea)
  library(pROC)
  library(ggplot2)
  library(ggpubr)
  library(batchtools)
  library(plyr)
  
  # source files
  source("batchtools_helper.R")
  
  pws <- c("shared", "independent")
  coverage <- c("real", 0.3, 1)
  rho <- c(0.2, 0.3, 0.8)
  
  batch.table <- expand.grid(pws,
                             coverage,
                             rho,
                             #type,
                             stringsAsFactors = F)
  batch.table$label <- paste0("cor",
                              "_rho", batch.table$Var3,
                              "_cov", batch.table$Var2,
                              "_", batch.table$Var1)
  batch.table$job <- paste0("cor",
                            "_rho", batch.table$Var3,
                            "_cov", batch.table$Var2,
                            "_", batch.table$Var1)#,
  #"_", batch.table$Var4)
  names <- paste0("eval", 1:dim(batch.table)[1])
  
  
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
  
  label = batch.table$label[kkk]
  #type = batch.table$Var4[kkk]
  outdir = "../results/current/thesis"
  iter = 1:100
  
  types <- c("GSEA_P", "GSEA_T",
             "MGSA_P", "MGSA_T",
             "MONA_1D", "MONA_1Dpm",
             "MONA_2D", "MONA_2Dpm")
  
  calc_perf <- function(sim.act, sim.inact, act, inact, roc=NULL){
    
    TP <- sum(sim.act %in% act)
    TN <- sum(sim.inact %in% inact)
    FP <- sum(sim.act %in% inact)
    FN <- sum(sim.inact %in% act)
    
    EXACT <- isTRUE(all.equal(sim.act, act))
    ACC <- (TP + TN)/(TP + FN + FP + TN)
    F1 <- (2*TP)/(2*TP + FP + FN)
    FPR <- (FP)/(FP + TN)
    FNR <- (FN)/(FN + TP)
    FDR <- ifelse(TP+FP==0, 0, (FP)/(TP + FP))
    SENS <- (TP)/(TP + FN)
    SPEC <- (TN)/(TN + FP)
    PREC<- ifelse(TP+FP==0, 0, (TP)/(TP + FP))
    REC <- (TP)/(TP + FN)
    
    AUC <- ifelse(is.null(roc), NA, as.numeric(auc(roc)))
    
    res <- as.numeric(c(TP, TN, FP, FN,
                        EXACT, ACC, F1,
                        FPR, FNR, FDR,
                        SENS, SPEC, PREC, REC,
                        AUC))
    names(res) <- c("TP", "TN", "FP", "FN",
                    "EXACT", "ACC", "F1",
                    "FPR", "FNR", "FDR",
                    "SENS", "SPEC", "PREC", "REC",
                    "AUC")
    return(res)
  }
  
  if (!dir.exists(paste0(outdir, "/", label, "/eval/"))) {
    dir.create(paste0(outdir, "/", label, "/eval/"), recursive = T)
  }
  
  for(k in iter){
    #k = 1
    sim_all <- readRDS(file=paste0(outdir, "/", label, "/sum_stats/",
                                   "sim_data_all_", k, ".RDS"))
    meta.data <- sim_all[[paste0(error.rates[1, "alpha"], "_", error.rates[1, "beta"])]]$meta.data
    
    gmt.file <- "../data/current/gmt_files/c2.cp.kegg.v7.2.symbols.gmt"
    gmt <- gmtPathways(gmt.file)
    active.genes <- meta.data$active.genes$gene
    
    ov <- data.frame(matrix(NA, dimnames = list(names(gmt),
                                                c("pathway", "overlap")),
                            nrow = length(gmt), ncol = 2),
                     stringsAsFactors = F)
    for (pw in 1:dim(ov)[1]){
      ov[pw, "pathway"] <- names(gmt)[pw]
      ov[pw, "overlap"] <- as.numeric(sum(active.genes %in% gmt[[pw]])/length(gmt[[pw]]))
    }
    active.terms.ov25 <- ov[ov$overlap>=0.25, "pathway"]
    
    eval_res <- data.frame(matrix(NA, 0, 21),
                           stringsAsFactors = F)
    colnames(eval_res) <- c("scenario", "type", "eval", "alpha", "beta", "N",
                            "TP", "TN", "FP", "FN",
                            "EXACT", "ACC", "F1",
                            "FPR", "FNR", "FDR",
                            "SENS", "SPEC", "PREC", "REC",
                            "AUC")
    eval_act_gt <- list()
    eval_act_sim <- list()
    roc_obj <- list()
    r <- 1
    
    for(type in types){
      # type <- types[1]
      
      if(type %in% c("GSEA_P", "GSEA_T")){
        res <- readRDS(paste0(outdir, "/", label, "/GSEA/", type, "_res/", type, "_res_", k, ".RDS"))
        res$pathway <- as.character(res$pathway)
        res <- res[!((is.na(res$pval_spec1) & is.na(res$pval_spec2)) | (is.na(res$pathway))), ]
      }
      
      if(type %in% c("MGSA_P", "MGSA_T")){
        res <- readRDS(paste0(outdir, "/", label, "/MGSA/", type, "_res/", type, "_res_", k, ".RDS"))
        res$pathway <- as.character(res$pathway)
        res <- res[!((is.na(res$posterior_spec1) & is.na(res$posterior_spec2)) | (is.na(res$pathway))), ]
      }
      
      if(type %in% c("MONA_1D", "MONA_1Dpm", "MONA_2D", "MONA_2Dpm")){
        res <- readRDS(paste0(outdir, "/", label, "/MONA/", type, "_res/", type, "_res_", k, ".RDS"))
        res$pathway <- as.character(res$pathway)
        if(type %in% c("MONA_1D", "MONA_1Dpm")){
          res <- res[!((is.na(res$posterior_spec1) & is.na(res$posterior_spec2)) | (is.na(res$pathway))), ]
        }
        if(type %in% c("MONA_2D", "MONA_2Dpm")){
          res <- res[!(is.na(res$posterior_multi) | is.na(res$pathway)), ]
        }
      }
      
      for (l in 1:dim(error.rates)[1]){
        # l=1
        alpha <- error.rates[l, "alpha"]
        beta <- error.rates[l, "beta"]
        res2 <- res[res$alpha==alpha & res$beta==beta & res$N==k, ]
        rownames(res2) <- res2$pathway
        
        evals = c("both",
                  "both.spec1", "both.spec2",
                  "both.AND", #"both.AND.sign",
                  "both.OR", #"both.OR.sign",
                  "spec1.spec1", "spec1.both",
                  "spec2.spec2", "spec2.both",
                  "both.strict",
                  "both.spec1.strict", "both.spec2.strict",
                  "both.AND.strict",
                  "both.OR.strict",
                  "spec1.spec1.strict", "spec1.both.strict",
                  "spec2.spec2.strict", "spec2.both.strict",
                  "both.AND.ov25", "both.OR.ov25",
                  "both.AND.top", "both.OR.top")
        
        for (eval in evals){
          # eval = evals[1]
          eval_temp <- data.frame("scenario" = label,
                                  "type" = type,
                                  "eval" = eval,
                                  "alpha" = alpha,
                                  "beta" = beta,
                                  "N" = k,
                                  stringsAsFactors = F)
          
          if(eval=="both.AND" & type %in% c("GSEA_P", "GSEA_T",
                                            "MGSA_P", "MGSA_T",
                                            "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim[which(res2$ES_spec1>0 & res2$ES_spec2>0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1>0 & res2$ES_spec2>0)], "_up")
                res2$act_sim[which(res2$ES_spec1<0 & res2$ES_spec2<0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1<0 & res2$ES_spec2<0)], "_down")
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.05 & res2$padj_spec2<0.05),
                                   "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)],
                                           res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)],
                                                   res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species=c(rep("spec1", sum(!is.na(res2$pval_spec1))),
                                              rep("spec2", sum(!is.na(res2$pval_spec2)))),
                                 stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(gsub("_up", "",
                                   gsub("_down", "",
                                        res2[which(res2$posterior_spec1>0.5 & res2$posterior_spec2>0.5),
                                             "pathway"])))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)],
                                                        res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)],
                                           res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species=c(rep("spec1", sum(!is.na(res2$posterior_spec1))),
                                              rep("spec2", sum(!is.na(res2$posterior_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc[, c("pathway", "species")]), ]
              res.roc <- res.roc[order(res.roc$score), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec1>0.5 &
                                           res2$posterior_spec2>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="both.OR" & type %in% c("GSEA_P", "GSEA_T",
                                           "MGSA_P", "MGSA_T",
                                           "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim[which(res2$ES_spec1>0 & res2$ES_spec2>0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1>0 & res2$ES_spec2>0)], "_up")
                res2$act_sim[which(res2$ES_spec1<0 & res2$ES_spec2<0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1<0 & res2$ES_spec2<0)], "_down")
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.05 | res2$padj_spec2<0.05),
                                   "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)],
                                              res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)],
                                                   res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species=c(rep("spec1", sum(!is.na(res2$pval_spec1))),
                                              rep("spec2", sum(!is.na(res2$pval_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec1>0.5 | res2$posterior_spec2>0.5),
                                                    "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)],
                                                        res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)],
                                            res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species=c(rep("spec1", sum(!is.na(res2$posterior_spec1))),
                                              rep("spec2", sum(!is.na(res2$posterior_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec1>0.5 |
                                           res2$posterior_spec2>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="both" & type %in% c("MONA_2D", "MONA_2Dpm")){
            if(type %in% c("MONA_2D", "MONA_2Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_multi>0.5), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      res2$pathway[!is.na(res2$posterior_multi)])),
                                    score=res2$posterior_multi[!is.na(res2$posterior_multi)],
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_multi>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="both.spec1" & type %in% c("MONA_2D", "MONA_2Dpm")){
            if(type %in% c("MONA_2D", "MONA_2Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_multi>0.5), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      res2$pathway[!is.na(res2$posterior_multi)])),
                                    score=res2$posterior_multi[!is.na(res2$posterior_multi)],
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_multi>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$act.terms$terms.spec1)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="both.spec2" & type %in% c("MONA_2D", "MONA_2Dpm")){
            if(type %in% c("MONA_2D", "MONA_2Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_multi>0.5), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      res2$pathway[!is.na(res2$posterior_multi)])),
                                    score=res2$posterior_multi[!is.na(res2$posterior_multi)],
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_multi>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$act.terms$terms.spec2)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="spec1.spec1" & type %in% c("GSEA_P", "GSEA_T",
                                               "MGSA_P", "MGSA_T",
                                               "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2[!is.na(res2$padj_spec1), "pathway"])
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim <- paste0(res2$pathway,
                                       ifelse(res2$ES_spec1>0,
                                              "_up", "_down"))
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.05), "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)])),
                                    species="spec1",
                                    stringsAsFactors = F)
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "",
                                   gsub("_down", "",
                                        res2[!is.na(res2$posterior_spec1), "pathway"])))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec1>0.5), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)]),
                                    species="spec1",
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec1>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$act.terms$terms.spec1)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="spec1.both" & type %in% c("GSEA_P", "GSEA_T",
                                              "MGSA_P", "MGSA_T",
                                              "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim <- paste0(res2$pathway,
                                       ifelse(res2$ES_spec1>0,
                                              "_up", "_down"))
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.05), "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)])),
                                    species="spec1",
                                    stringsAsFactors = F)
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "",
                                   gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec1>0.5), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)]),
                                    species="spec1",
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec1>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="spec2.spec2" & type %in% c("GSEA_P", "GSEA_T",
                                               "MGSA_P", "MGSA_T",
                                               "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2[!is.na(res2$padj_spec2), "pathway"])
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim <- paste0(res2$pathway,
                                       ifelse(res2$ES_spec2>0,
                                              "_up", "_down"))
              }
              sim.act <- sort(res2[which(res2$padj_spec2<0.05), "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species="spec2",
                                    stringsAsFactors = F)
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "",
                                   gsub("_down", "",
                                        res2[!is.na(res2$posterior_spec2), "pathway"])))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec2>0.5), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species="spec2",
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec2>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$act.terms$terms.spec2)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="spec2.both" & type %in% c("GSEA_P", "GSEA_T",
                                              "MGSA_P", "MGSA_T",
                                              "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim <- paste0(res2$pathway,
                                       ifelse(res2$ES_spec2>0,
                                              "_up", "_down"))
              }
              sim.act <- sort(res2[which(res2$padj_spec2<0.05), "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species="spec2",
                                    stringsAsFactors = F)
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "",
                                   gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec2>0.5), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species="spec2",
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec2>0.5), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="both.AND.strict" & type %in% c("GSEA_P", "GSEA_T",
                                                   "MGSA_P", "MGSA_T",
                                                   "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim[which(res2$ES_spec1>0 & res2$ES_spec2>0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1>0 & res2$ES_spec2>0)], "_up")
                res2$act_sim[which(res2$ES_spec1<0 & res2$ES_spec2<0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1<0 & res2$ES_spec2<0)], "_down")
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.01 & res2$padj_spec2<0.01),
                                   "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)],
                                              res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)],
                                                   res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species=c(rep("spec1", sum(!is.na(res2$pval_spec1))),
                                              rep("spec2", sum(!is.na(res2$pval_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(gsub("_up", "",
                                   gsub("_down", "",
                                        res2[which(res2$posterior_spec1>0.6 & res2$posterior_spec2>0.6),
                                             "pathway"])))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)],
                                                        res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)],
                                            res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species=c(rep("spec1", sum(!is.na(res2$posterior_spec1))),
                                              rep("spec2", sum(!is.na(res2$posterior_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc[, c("pathway", "species")]), ]
              res.roc <- res.roc[order(res.roc$score), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec1>0.6 &
                                           res2$posterior_spec2>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="both.OR.strict" & type %in% c("GSEA_P", "GSEA_T",
                                                  "MGSA_P", "MGSA_T",
                                                  "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim[which(res2$ES_spec1>0 & res2$ES_spec2>0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1>0 & res2$ES_spec2>0)], "_up")
                res2$act_sim[which(res2$ES_spec1<0 & res2$ES_spec2<0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1<0 & res2$ES_spec2<0)], "_down")
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.01 | res2$padj_spec2<0.01),
                                   "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)],
                                              res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)],
                                                   res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species=c(rep("spec1", sum(!is.na(res2$pval_spec1))),
                                              rep("spec2", sum(!is.na(res2$pval_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec1>0.6 | res2$posterior_spec2>0.6),
                                                    "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)],
                                                        res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)],
                                            res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species=c(rep("spec1", sum(!is.na(res2$posterior_spec1))),
                                              rep("spec2", sum(!is.na(res2$posterior_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec1>0.6 |
                                           res2$posterior_spec2>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="both.strict" & type %in% c("MONA_2D", "MONA_2Dpm")){
            if(type %in% c("MONA_2D", "MONA_2Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_multi>0.6), "pathway"]))))
              act_sim <- sort(res2[which(res2$posterior_multi>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            res.roc <- data.frame(pathway=gsub("_up", "",
                                               gsub("_down", "",
                                                    res2$pathway[!is.na(res2$posterior_multi)])),
                                  score=c(res2$posterior_multi[!is.na(res2$posterior_multi)]),
                                  stringsAsFactors = F)
            res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
            res.roc <- res.roc[!duplicated(res.roc$pathway), ]
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="both.spec1.strict" & type %in% c("MONA_2D", "MONA_2Dpm")){
            if(type %in% c("MONA_2D", "MONA_2Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_multi>0.6), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      res2$pathway[!is.na(res2$posterior_multi)])),
                                    score=res2$posterior_multi[!is.na(res2$posterior_multi)],
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_multi>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$act.terms$terms.spec1)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="both.spec2.strict" & type %in% c("MONA_2D", "MONA_2Dpm")){
            if(type %in% c("MONA_2D", "MONA_2Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_multi>0.6), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      res2$pathway[!is.na(res2$posterior_multi)])),
                                    score=res2$posterior_multi[!is.na(res2$posterior_multi)],
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_multi>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$act.terms$terms.spec2)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="spec1.spec1.strict" & type %in% c("GSEA_P", "GSEA_T",
                                                      "MGSA_P", "MGSA_T",
                                                      "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2[!is.na(res2$padj_spec1), "pathway"])
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim <- paste0(res2$pathway,
                                       ifelse(res2$ES_spec1>0,
                                              "_up", "_down"))
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.01), "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)])),
                                    species="spec1",
                                    stringsAsFactors = F)
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "",
                                   gsub("_down", "",
                                        res2[!is.na(res2$posterior_spec1), "pathway"])))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec1>0.6), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)]),
                                    species="spec1",
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec1>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$act.terms$terms.spec1)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="spec1.both.strict" & type %in% c("GSEA_P", "GSEA_T",
                                                     "MGSA_P", "MGSA_T",
                                                     "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim <- paste0(res2$pathway,
                                       ifelse(res2$ES_spec1>0,
                                              "_up", "_down"))
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.01), "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)])),
                                    species="spec1",
                                    stringsAsFactors = F)
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "",
                                   gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec1>0.6), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)]),
                                    species="spec1",
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec1>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="spec2.spec2.strict" & type %in% c("GSEA_P", "GSEA_T",
                                                      "MGSA_P", "MGSA_T",
                                                      "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2[!is.na(res2$padj_spec2), "pathway"])
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim <- paste0(res2$pathway,
                                       ifelse(res2$ES_spec2>0,
                                              "_up", "_down"))
              }
              sim.act <- sort(res2[which(res2$padj_spec2<0.01), "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species="spec2",
                                    stringsAsFactors = F)
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "",
                                   gsub("_down", "",
                                        res2[!is.na(res2$posterior_spec2), "pathway"])))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec2>0.6), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species="spec2",
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec2>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$act.terms$terms.spec2)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="spec2.both.strict" & type %in% c("GSEA_P", "GSEA_T",
                                                     "MGSA_P", "MGSA_T",
                                                     "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim <- paste0(res2$pathway,
                                       ifelse(res2$ES_spec2>0,
                                              "_up", "_down"))
              }
              sim.act <- sort(res2[which(res2$padj_spec2<0.01), "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species="spec2",
                                    stringsAsFactors = F)
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "",
                                   gsub("_down", "", res2$pathway)))
              sim.act <- sort(unique(gsub("_up", "",
                                          gsub("_down", "",
                                               res2[which(res2$posterior_spec2>0.6), "pathway"]))))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species="spec2",
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[which(res2$posterior_spec2>0.6), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="both.AND.ov25" & type %in% c("GSEA_P", "GSEA_T")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim[which(res2$ES_spec1>0 & res2$ES_spec2>0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1>0 & res2$ES_spec2>0)], "_up")
                res2$act_sim[which(res2$ES_spec1<0 & res2$ES_spec2<0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1<0 & res2$ES_spec2<0)], "_down")
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.05 & res2$padj_spec2<0.05),
                                   "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)],
                                              res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)],
                                                   res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species=c(rep("spec1", sum(!is.na(res2$pval_spec1))),
                                              rep("spec2", sum(!is.na(res2$pval_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(active.terms.ov25)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act, "_",
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "up", "down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="both.OR.ov25" & type %in% c("GSEA_P", "GSEA_T")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim[which(res2$ES_spec1>0 & res2$ES_spec2>0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1>0 & res2$ES_spec2>0)], "_up")
                res2$act_sim[which(res2$ES_spec1<0 & res2$ES_spec2<0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1<0 & res2$ES_spec2<0)], "_down")
              }
              sim.act <- sort(res2[which(res2$padj_spec1<0.05 | res2$padj_spec2<0.05),
                                   "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)],
                                              res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)],
                                                   res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species=c(rep("spec1", sum(!is.na(res2$pval_spec1))),
                                              rep("spec2", sum(!is.na(res2$pval_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(active.terms.ov25)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act, "_",
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "up", "down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          if(eval=="both.AND.top" & type %in% c("GSEA_P", "GSEA_T",
                                                "MGSA_P", "MGSA_T",
                                                "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim[which(res2$ES_spec1>0 & res2$ES_spec2>0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1>0 & res2$ES_spec2>0)], "_up")
                res2$act_sim[which(res2$ES_spec1<0 & res2$ES_spec2<0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1<0 & res2$ES_spec2<0)], "_down")
              }
              sim.act <- sort(res2[which(res2$pval_spec1 <= sort(res2$pval_spec1)[6] &
                                           res2$pval_spec2 <= sort(res2$pval_spec2)[6]),
                                   "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)],
                                              res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)],
                                                   res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species=c(rep("spec1", sum(!is.na(res2$pval_spec1))),
                                              rep("spec2", sum(!is.na(res2$pval_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              res3 <- res2[which(((res2$posterior_spec1 >= min(sort(res2$posterior_spec1,
                                                                    decreasing = T)[1:6])) &
                                    (res2$posterior_spec1>0)) &
                                   ((res2$posterior_spec2 >= min(sort(res2$posterior_spec2,
                                                                      decreasing = T)[1:6])) &
                                      (res2$posterior_spec2>0))), ]

              sim.act <- sort(gsub("_up", "",
                                   gsub("_down", "",
                                        res3[which((res3$posterior_spec1 > 0) &
                                                     (res3$posterior_spec2 > 0)), "pathway"])))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)],
                                                        res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)],
                                            res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species=c(rep("spec1", sum(!is.na(res2$posterior_spec1))),
                                              rep("spec2", sum(!is.na(res2$posterior_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc[, c("pathway", "species")]), ]
              res.roc <- res.roc[order(res.roc$score), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res3[which((res3$posterior_spec1 > 0) &
                                           (res3$posterior_spec2 > 0)), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          if(eval=="both.OR.top" & type %in% c("GSEA_P", "GSEA_T",
                                               "MGSA_P", "MGSA_T",
                                               "MONA_1D", "MONA_1Dpm")){
            if(type %in% c("GSEA_P", "GSEA_T")){
              terms <- unique(res2$pathway)
              res2$act_sim <- paste0(res2$pathway, "_NA")
              if(type=="GSEA_T"){
                res2$act_sim[which(res2$ES_spec1>0 & res2$ES_spec2>0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1>0 & res2$ES_spec2>0)], "_up")
                res2$act_sim[which(res2$ES_spec1<0 & res2$ES_spec2<0)] <-
                  paste0(res2$pathway[which(res2$ES_spec1<0 & res2$ES_spec2<0)], "_down")
              }
              sim.act <- sort(res2[which(res2$pval_spec1 <= sort(res2$pval_spec1)[6] |
                                           res2$pval_spec2 <= sort(res2$pval_spec2)[6]),
                                   "pathway"])
              res.roc <- data.frame(pathway=c(res2$pathway[!is.na(res2$pval_spec1)],
                                              res2$pathway[!is.na(res2$pval_spec2)]),
                                    score=-log10(c(res2$pval_spec1[!is.na(res2$pval_spec1)],
                                                   res2$pval_spec2[!is.na(res2$pval_spec2)])),
                                    species=c(rep("spec1", sum(!is.na(res2$pval_spec1))),
                                              rep("spec2", sum(!is.na(res2$pval_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res2[sim.act, "act_sim"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            if(type %in% c("MGSA_P", "MGSA_T", "MONA_1D", "MONA_1Dpm")){
              terms <- unique(gsub("_up", "", gsub("_down", "", res2$pathway)))
              res3 <- res2[which(((res2$posterior_spec1 >= min(sort(res2$posterior_spec1,
                                                                    decreasing = T)[1:6])) &
                                    (res2$posterior_spec1>0)) |
                                   ((res2$posterior_spec2 >= min(sort(res2$posterior_spec2,
                                                                      decreasing = T)[1:6])) &
                                      (res2$posterior_spec2>0))), ]
              sim.act <- sort(gsub("_up", "",
                                   gsub("_down", "",
                                        res3[which((res3$posterior_spec1 > 0) |
                                                     (res3$posterior_spec2 > 0)), "pathway"])))
              res.roc <- data.frame(pathway=gsub("_up", "",
                                                 gsub("_down", "",
                                                      c(res2$pathway[!is.na(res2$posterior_spec1)],
                                                        res2$pathway[!is.na(res2$posterior_spec2)]))),
                                    score=c(res2$posterior_spec1[!is.na(res2$posterior_spec1)],
                                            res2$posterior_spec2[!is.na(res2$posterior_spec2)]),
                                    species=c(rep("spec1", sum(!is.na(res2$posterior_spec1))),
                                              rep("spec2", sum(!is.na(res2$posterior_spec2)))),
                                    stringsAsFactors = F)
              res.roc <- res.roc[order(res.roc$score, decreasing = T), ]
              res.roc <- res.roc[!duplicated(res.roc$pathway), ]
              act_sim <- sort(res3[which((res3$posterior_spec1 > 0) |
                                           (res3$posterior_spec2 > 0)), "pathway"])
              sim.inact <- terms[!(terms %in% sim.act)]
            }
            act <- sort(meta.data$active.terms$term)
            inact <- terms[!(terms %in% act)]
            
            roc <- NULL
            if(sum(res.roc$score==Inf)>0){
              res.roc$score[res.roc$score==Inf] <- max(res.roc$score[res.roc$score!=Inf]) + 1
            }
            roc <- roc(as.numeric(res.roc$pathway %in% act), res.roc$score,
                       direction = "<", levels = c("0", "1"))
            
            eval_temp[1, c("TP", "TN", "FP", "FN",
                           "EXACT", "ACC", "F1",
                           "FPR", "FNR", "FDR", "SENS", "SPEC", "PREC", "REC", "AUC")] <-
              calc_perf(sim.act, sim.inact, act, inact, roc)
            
            eval_res <- rbind(eval_res, eval_temp)
            eval_act_gt[[r]] <- list(paste0(act,
                                            ifelse(meta.data$active.terms[act, "sign"]>0,
                                                   "_up", "_down")))
            eval_act_sim[[r]] <- list(as.character(act_sim))
            roc_obj[[r]] <- roc
            r <- r+1
          }
          
          # if(eval=="spec1.spec1.sign"){
          #   sim.act <- res2[which(res2$padj_spec1<0.05 &
          #                           sign(res2$ES_spec1)==sign(meta.data$active.terms[res2$pathway, "sign"])),
          #                   "pathway"]
          #   sim.inact <- res2[which(!(res2$padj_spec1<0.05)),
          #                     "pathway"]
          #   act <- meta.data$act.terms$terms.spec1
          #   terms <- c(sim.act, sim.inact)
          #   inact <- terms[!(terms %in% act)]
          #   eval_temp <- cbind(calc_perf(sim.act, sim.inact, act, inact))
          #   eval_res <- rbind(eval_res, eval_temp)
          # }
        }
      }
    }
    
    result <- data.frame(eval_res,
                         I(eval_act_gt), I(eval_act_sim), I(roc_obj),
                         stringsAsFactors = F)
    
    saveRDS(result,
            file=paste0(outdir, "/", label, "/eval/",
                        "evaluation_", k, ".RDS"))
    
    print(paste0("Evaluation done for ", label, " (counter=", kkk, "), iteration ", k, " (of ", length(iter), ")"))
  }
}

pws <- c("shared", "independent")
coverage <- c("real", 0.3, 1)
rho <- c(0.2, 0.3, 0.8)

batch.table <- expand.grid(pws,
                           coverage,
                           rho,
                           #type,
                           stringsAsFactors = F)
batch.table$label <- paste0("cor",
                            "_rho", batch.table$Var3,
                            "_cov", batch.table$Var2,
                            "_", batch.table$Var1)
batch.table$job <- paste0("cor",
                          "_rho", batch.table$Var3,
                          "_cov", batch.table$Var2,
                          "_", batch.table$Var1)

names <- paste0("eval", 1:dim(batch.table)[1])

run.batchtools(run_eval, 1:dim(batch.table)[1], more.args=list(),
               names, "tmp_dir_eval_thesis",
               clean.up=F,
               resources=list(partition="cpu_p",
                              exclude="cpusrv07,cpusrv26,cpusrv28,ibis216-010-035,ibis216-010-037,ibis216-010-051,ibis216-010-064",
                              memory="20G",
                              ncpus = 1,
                              measure.memory = TRUE,
                              walltime="24:00:00"))



setwd("~/work/multi_omics_enrich/scripts/")
library(batchtools)
loadRegistry(file.dir="tmp_dir_eval_thesis", writeable = F)
getStatus
getJobTable()

if(length(c(findExpired()$job.id, findErrors()$job.id)>0)){
  loadRegistry(file.dir="tmp_dir_eval_thesis", writeable = T)
  getStatus()
  del_ids <- NULL
  del_ids <- c(findExpired()$job.id, findErrors()$job.id)
  if(!is.null(del_ids)){
    submitJobs(del_ids,
               resources=list(partition="cpu_p",
                              memory="20G",
                              ncpus = 1,
                              measure.memory = TRUE,
                              walltime="24:00:00"))
  }
  getStatus()
}

q(save = "no")
