#!/usr/bin/env Rscript
library(fgsea)
library(MASS)
library(pROC)
# library(svMisc)
### OptParse ########
library(optparse)
option_list <- list( 
  make_option(c( "--out"), type="character", default="eval.png", 
              help="Outputpath of image. (Default: eval.png"), 
  make_option(c("--pred"), type="character", default="pred.rds",
              help="Prediction file."),
  make_option(c("--gtruth"), type="character", default="gtruth.rds", 
              help="Simulation file."),
  make_option(c("--sign"), type="character", default="yes", 
              help="Include downregulated pathways? (Default: yes)")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


out     <- opt$out
sign <- if(opt$sign=="yes")
pred <- readRDS(opt$pred)
gtruth <- readRDS(opt$gtruth)
#####################
### Functions #######
eval_sim_single <- function(pred, gtruth, sign=F, out){

  
  if(sign){
    gtruth$active.terms <- .rename_terms(gtruth$active.terms) 
    gtruth$act.terms <- .rename_terms2(gtruth$act.terms)
  }
  
  predPWs <- as.character(pred$res[pred$res$posterior>0.6,1])
  truePWs.mRNA <- gtruth$act.terms$terms.trans
  truePWs.prot <- gtruth$act.terms$terms.prot
  
  # df.mRNA <- .evalDF(predPWs, truePWs.mRNA)
  # df.prot <- .evalDF(predPWs, truePWs.prot)
  # barplot(as.matrix(df.mRNA))
  
  
  
  
  TP.mRNA <- sum(as.numeric(predPWs %in% truePWs.mRNA))
  TP.prot <- sum(as.numeric(predPWs %in% truePWs.prot))
  FP.mRNA <- sum(as.numeric(!predPWs %in% truePWs.mRNA))
  FP.prot <- sum(as.numeric(!predPWs %in% truePWs.prot))
  FN.mRNA <- sum(as.numeric(!truePWs.mRNA %in% predPWs))
  FN.prot <- sum(as.numeric(!truePWs.prot %in% predPWs))
  hits <- data.frame(mRNA = c(TP.mRNA, FP.mRNA, FN.mRNA), prot = c(TP.prot, FP.prot, FN.prot))
  rownames(hits) <- c("TP", "FP", "FN")
  png(out)
  barplot(as.matrix(hits), col = c("green", "red", "orange"), main = "Evaluation", sub = paste0("Predicted ",TP.mRNA,"/",length(truePWs.mRNA)," in mRNA and ",TP.prot,"/",length(truePWs.prot)," in proteins correctly."))
  legend("topleft",
         c("TP","FP", "FN"),
         fill = c("green", "red", "orange")
  )
  dev.off()
}
.evalDF <- function(truePWs, predPWs){
  data.frame(pathways = unique(predPWs, truePWs), predicted = as.numeric(unique(predPWs, truePWs) %in% predPWs), true = as.numeric(unique(predPWs, truePWs) %in% truePWs))
}
.rename_terms <- function(activePWsign){
  for (row in c(1:nrow(activePWsign))) {
    if(activePWsign[row,2]<0)  activePWsign[row,1] = paste0(activePWsign[row,1],"_down")
    if(activePWsign[row,2]>0) activePWsign[row,1] = paste0(activePWsign[row,1],"_up")
  }
  return(activePWsign)
}
.rename_terms2 <- function(activePWsign){
  for (row in c(1:nrow(activePWsign))) {
    if(activePWsign[row,2]<0)  activePWsign[row,1] = paste0(activePWsign[row,1],"_down")
    if(activePWsign[row,2]>0) activePWsign[row,1] = paste0(activePWsign[row,1],"_up")
    if(activePWsign[row,4]<0)  activePWsign[row,3] = paste0(activePWsign[row,3],"_down")
    if(activePWsign[row,4]>0) activePWsign[row,3] = paste0(activePWsign[row,3],"_up")
  }
  return(activePWsign)
}
#####################
### Main ############

eval_sim_single(pred, gtruth, sign, out)
#####################


