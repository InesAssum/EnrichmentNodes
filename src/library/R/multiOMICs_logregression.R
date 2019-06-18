#!/usr/bin/env Rscript
###############################
library(fgsea)
library(optparse)
library(pROC)
######### args ################
option_list = list(
  make_option("--pheno", type="character", default=NULL, help="Phenotype file (rows = samples, cols = covs)"),  
  make_option("--species1", type="character", default=NULL,help="Species 1 data (rows = samples, cols = genes)"),
  make_option("--species2", type="character", default=NULL,help="Species 2 data (rows = samples, cols = genes)"),
  make_option("--trait", type="character", default=NULL,help="Trait of interest (Must be column in phenotype file)"),
  make_option("--out", type="character", default="out.RDS", help="output file name [default= out.RDS]")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
pheno   <- opt$pheno
spec1   <- opt$species1
spec2   <- opt$species2
out     <- opt$out
trait   <- opt$trait
###############################
######### functions ###########
mcFaddens <- function(mod, nullmod){
  1-logLik(mod)/logLik(nullmod)
}

my_logreg <- function(pheno, spec1, spec2, traitofinterest){
  print("Reading data...")
  phenos   <- readRDS(pheno)
  expr   <- readRDS(spec1)
  prot   <- readRDS(spec2)
  covariates <- colnames(phenos)
  trait <- traitofinterest
  if(!(trait %in% covariates)) stop("Could not find trait of interest in Phenotypes file")
  
  # filter
  expr        <- expr[rownames(phenos),]
  prot        <- prot[rownames(phenos),]
  covariates  <- covariates[covariates != trait]
  print("Initialising...")
  # null test:
  df <- data.frame(phenos)
  df <- df[complete.cases(df),] 
  nullmod <- glm(data=df, family = binomial, formula = as.formula(paste0(trait," ~ ",paste(covariates,collapse = " + "))), na.action = na.exclude)
  auc_nullm = auc(roc(nullmod$y, nullmod$fitted.values)) # AUC
  
  # init
  genes <- unique(c(colnames(expr),colnames(prot)))
  result <- matrix(nrow = length(genes), ncol = 10)
  rownames(result) <- genes
  colnames(result) <- c('deviance', 'AIC', 'diff(AIC)','AUC','diff(AUC)', 'beta.species1', 'p-val.species1','beta.species2', 'p-val.species2', 'pseudo.R2')
  print("Done")
  print("Starting logistic regression...")
  for (i in 1:length(genes)) {
    gene <- genes[i]
    in_expr <- gene %in% colnames(expr)
    in_prot <- gene %in% colnames(prot)
    # case: only trans
    if(in_expr && !(in_prot)){
      df <- data.frame(trans=expr[,gene], phenos)
      df <- df[complete.cases(df),]
      formel1 <- paste0("overall_AF ~ ",paste(c("trans",covariates),collapse = " + "))
      glm1 <- glm(data=df, family = binomial, formula = as.formula(formel1), na.action = na.exclude)
      auc1 = auc(roc(glm1$y, glm1$fitted.values)) # AUC
      result[i,"deviance"]    <-glm1$deviance    # deviance
      result[i,"AIC"]         <-AIC(glm1)    # AIC
      result[i,"diff(AIC)"]   <-AIC(nullmod)-AIC(glm1)    # diff(AIC)
      result[i,"AUC"]         <-  auc1              # AUC
      result[i,"diff(AUC)"]   <-  auc_nullm-auc1    # diff(AUC)
      result[i,"beta.species1"]  <-coefficients(summary(glm1))[2,1]    # Estimate trans
      result[i,"p-val.species1"] <-coefficients(summary(glm1))[2,4]    # p-val trans
      result[i,"beta.species2"]   <-NA   # Estimate prot
      result[i,"p-val.species2"]  <-NA    # p-val prot
      result[i,"pseudo.R2"]   <-mcFaddens(glm1,nullmod)    # McFadden's R squared  
    }
    # case: only prot
    if(!(in_expr) && in_prot){
      df <- data.frame(prot=prot[,gene], phenos)
      df <- df[complete.cases(df),]
      formel1 <- paste0("overall_AF ~ ",paste(c("prot",covariates),collapse = " + "))
      glm1 <- glm(data=df, family = binomial, formula = as.formula(formel1), na.action = na.exclude)
      auc1 = auc(roc(glm1$y, glm1$fitted.values)) # AUC
      result[i,"deviance"]    <-glm1$deviance    # deviance
      result[i,"AIC"]         <-AIC(glm1)    # AIC
      result[i,"diff(AIC)"]   <-AIC(nullmod)-AIC(glm1)    # diff(AIC)
      result[i,"AUC"]         <-  auc1              # AUC
      result[i,"diff(AUC)"]   <-  auc_nullm-auc1    # diff(AUC)
      result[i,"beta.species1"]  <-NA    # Estimate trans
      result[i,"p-val.species1"] <-NA    # p-val trans
      result[i,"beta.species2"]   <-coefficients(summary(glm1))[2,1]    # Estimate prot
      result[i,"p-val.species2"]  <-coefficients(summary(glm1))[2,4]    # p-val prot
      result[i,"pseudo.R2"]   <-mcFaddens(glm1,nullmod)    # McFadden's R squared  
    }
    # case: no missing gene
    if(in_expr && in_prot){
      df <- data.frame(trans=expr[,gene], prot=prot[,gene], phenos)
      df <- df[complete.cases(df),]
      formel1 <- paste0("overall_AF ~ trans + prot + ",paste(covariates,collapse = " + "))
      glm1 <- glm(data=df, family = binomial, formula = as.formula(formel1), na.action = na.exclude)
      auc1 = auc(roc(glm1$y, glm1$fitted.values)) # AUC
      result[i,"deviance"]    <-glm1$deviance    # deviance
      result[i,"AIC"]         <-AIC(glm1)    # AIC
      result[i,"diff(AIC)"]   <-AIC(glm1)-AIC(nullmod)    # diff(AIC)
      result[i,"AUC"]         <-  auc1              # AUC
      result[i,"diff(AUC)"]   <-  auc_nullm-auc1    # diff(AUC)
      result[i,"beta.species1"]  <-coefficients(summary(glm1))[1,1]    # Estimate trans
      result[i,"p-val.species1"] <-coefficients(summary(glm1))[1,4]    # p-val trans
      result[i,"beta.species2"]   <-coefficients(summary(glm1))[2,1]    # Estimate prot
      result[i,"p-val.species2"]  <-coefficients(summary(glm1))[2,4]    # p-val prot
      result[i,"pseudo.R2"]   <-mcFaddens(glm1,nullmod)    # McFadden's R squared  
    }
  }
  print("Sucessful.")
  return(result)
}
###############################
########### MAIN ##############
saveRDS(my_logreg(pheno, spec1, spec2, trait), out)

###############################
########## DEBUGGING ##########
pheno           <- "Documents/AFHRI_B/AFHRI_B/phenotypes.RDS"
# spec1           <- "Documents/AFHRI_B/AFHRI_B/expr_example_short.RDS"
# spec2           <- "Documents/AFHRI_B/AFHRI_B/prot_example.RDS"
# out             <- "Documents/AFHRI_B/AFHRI_B/test_out.RDS"
# trait <- "overall_AF"
# result1 <- readRDS("Documents/AFHRI_B/AFHRI_B/test_out2.RDS")
###############################

