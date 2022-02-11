#' #############################################################################
# Simulate paired summary statistics for two omics  ############################
#' #############################################################################
#'
#'
#' Simulate paired summary statistics for two omics
#'
#' @param label label for the simulated data (string)
#' @param outdir path to out directory (string)
#' @param seed if NA, random seed chosen and documented
#'
#' @param iterations number of iterations
#' @param error.rates data.frame with fpr / fnr combinations
#'
#' @param gmt pathway definition (list)
#' @param gmt.path path to load gmt file (default NA)
#' @param coverage fraction of features missing in the second omic (numeric) or list with two vectors of characters
#' @param scenario string: "shared" or "independent"
#' @param nactive number of active pathways (both or each omic)
#' @param rho correlation
#' @param filter.miss remove pathways with less than x genes
#'
#' @author Ines Assum


simulate_summary_stats <- function(label, outdir,
                                   seed=NA,
                                   iterations=1:100, error.rates,
                                   gmt=NULL, gmt.path=NA,
                                   coverage=0.3, coverage.path=NA,
                                   scenario="shared", nactive = 6,
                                   rho=0.3, filter.miss=10,
                                   sd.bg=2,
                                   sd.signal="1.5-alpha"){
  
  # load libraries
  require("MASS")
  
  if(is.na(seed)){
    seed <- sample(0:9999, 1)
  }
  set.seed(seed)
  
  parameters <- data.frame(label=label,
                           outdir=outdir,
                           gmt=gmt.path,
                           coverage=coverage,
                           coverage.path=coverage.path,
                           scenario=scenario,
                           nactive=paste(unlist(nactive), collapse = ","),
                           rho=rho,
                           filter.miss=filter.miss,
                           sd.bg=sd.bg,
                           sd.signal=sd.signal,
                           seed=seed,
                           time.of.simulation=NA,
                           max.error=NA)
  print(parameters)
  
  # runtime tracking ----
  time <- system.time({
    
    # create out directories
    dir.create(paste0(outdir, "/", label, "/sum_stats/"),
               recursive = T, showWarnings = F)
    dir.create(paste0(outdir, "/", label, "/QC_stats/"),
               recursive = T, showWarnings = F)
    
    # prep pathway definition
    # filter pathways by measured transcripts / proteins
    if(is.null(gmt) & is.character(gmt.path)){
      require(fgsea)
      gmt <- gmtPathways(gmt.path)
    }else if(is.list(gmt)){
    }else{
      stop("Error: GMT parameter not appropriate")
    }
    hidden <- unique(unlist(gmt))
    Adj <- matrix(NA, dimnames = list(hidden, names(gmt)),
                  nrow = length(hidden), ncol = length(gmt))
    for (i in 1:dim(Adj)[2]){
      Adj[,i] <- as.numeric(hidden %in% gmt[[i]])
    }
    
    # filter pathways by measured transcripts / proteins
    if(is.character(coverage) & is.character(coverage.path)){
      coverage <- readRDS(coverage.path)
      spec1 <- as.character(intersect(hidden, c(coverage[["spec1"]], coverage[["spec2"]])))
      spec2 <- as.character(intersect(hidden, coverage[["spec2"]]))
    }else if(is.numeric(coverage)){
      spec1 <- hidden
      spec2 <- sample(hidden, coverage*length(hidden))
    }else{
      stop("Error: Coverage parameter not appropriate")
    }
    
    hidden2 <- unique(spec1, spec2)
    
    if(filter.miss){
      Adj.spec1 <- Adj[spec1, colnames(Adj)[colSums(Adj[spec1,])>filter.miss]]
      Adj.spec2 <- Adj[spec2, colnames(Adj)[colSums(Adj[spec2,])>filter.miss]]
      
      pw.spec1 <- colnames(Adj.spec1)
      pw.spec2 <- colnames(Adj.spec2)
    }
    
    # initialize data sets and QC stats
    qc.spec1 <- data.frame(matrix(NA,0,7))
    qc.spec2 <- data.frame(matrix(NA,0,7))
    
    #TODO: pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) ----
    
    # for each iteration ---------------------------------------------------------
    for (k in iterations){
      # k=1
      sim.data <- list()
      
      ### sample terms, genes and signs --------------------------------------
      # sample active terms and their sign
      if(scenario=="independent"){
        act.terms <- data.frame(terms.spec1=as.character(sample(pw.spec1,
                                                                nactive[["spec1"]])),
                                signs.spec1=NA,
                                terms.spec2=as.character(sample(pw.spec2,
                                                                nactive[["spec2"]])),
                                signs.spec2=NA,
                                stringsAsFactors = F)
        active.terms <- unique(c(act.terms$terms.spec1,
                                 act.terms$terms.spec2))
        active.terms <- data.frame(term=active.terms,
                                   sign=sample(c(-1,1), length(active.terms),replace=T),
                                   row.names = active.terms,
                                   stringsAsFactors = F)
        
        act.terms$signs.spec1 <- active.terms[act.terms$terms.spec1, "sign"]
        act.terms$signs.spec2 <- active.terms[act.terms$terms.spec2, "sign"]
        
        # print.data.frame(act.terms,
        #                  row.names = F)
        
      }else if(scenario=="shared"){
        active.terms <- as.character(sample(unique(#pw.spec1,
          pw.spec2),
          nactive[["shared"]]))
        act.terms <- data.frame(terms.spec1=active.terms,
                                signs.spec1=NA,
                                terms.spec2=active.terms,
                                signs.spec2=NA,
                                stringsAsFactors = F)
        
        active.terms <- data.frame(term=active.terms,
                                   sign=sample(c(-1,1), length(active.terms),replace=T),
                                   row.names = active.terms,
                                   stringsAsFactors = F)
        
        act.terms$signs.spec1 <- active.terms[act.terms$terms.spec1, "sign"]
        act.terms$signs.spec2 <- active.terms[act.terms$terms.spec2, "sign"]
        
        # print.data.frame(act.terms,
        #                  row.names = F)
      }
      
      # assign regulated genes and their sign
      # (depends on sign of all active terms it belongs to, since a gene can
      #  belong to multiple pathways)
      act.genes.spec1 <- data.frame(gene = rownames(Adj.spec1)[rowSums(Adj.spec1[spec1,
                                                                                 act.terms$terms.spec1])>0],
                                    sign = NA,
                                    stringsAsFactors = F)
      rownames(act.genes.spec1) <- act.genes.spec1$gene
      act.genes.spec2 <- data.frame(gene = rownames(Adj.spec2)[rowSums(Adj.spec2[spec2,
                                                                                 act.terms$terms.spec2])>0],
                                    sign = NA,
                                    stringsAsFactors = F)
      rownames(act.genes.spec2) <- act.genes.spec2$gene
      active.genes <- data.frame(gene = as.character(unique(c(act.genes.spec1$gene, act.genes.spec2$gene))),
                                 sign = NA,
                                 stringsAsFactors = F)
      rownames(active.genes) <- active.genes$gene
      for (m in 1:length(active.genes$gene)){
        active.genes$sign[m] <-
          sample(active.terms$sign[as.logical(Adj[active.genes$gene[m], active.terms$term])], 1)
      }
      rownames(active.genes) <- active.genes$gene
      
      act.genes.spec1$sign <- active.genes[act.genes.spec1$gene, "sign"]
      act.genes.spec2$sign <- active.genes[act.genes.spec2$gene, "sign"]
      
      print("Active pathways:")
      print.data.frame(active.terms)
      
      ## for each combination of alpha & beta ------------------------------------
      for (l in 1:dim(error.rates)[1]){
        
        # # set seed to make simulations reproducible
        # seed <- sample(1:1e6, 1)
        # set.seed(seed)
        # print(paste0("Seed set to ", seed, "."))
        
        # l=1
        alpha <- error.rates[l, "alpha"]
        beta <- error.rates[l, "beta"]
        
        # set sigma of normal distribution for background and signal based on alpha and beta
        if(sd.bg=="1+100*beta"){
          sd.bg <- 1+100*beta
        }else if(is.numeric(sd.bg)){
          sd.bg <- sd.bg
        }else if(is.character(sd.bg)){
          error("Your term for the background standard deviation is not yet supported.")
        }
        
        if(sd.signal=="(1-alpha)*beta"){
          sd.sig <- (1-alpha)*beta
        }else if(sd.signal=="(1-alpha)"){
          sd.sig <- (1-alpha)
        }else if(sd.signal=="(1-alpha)*100*beta"){
          sd.sig <- (1-alpha)*100*beta
        }else if(sd.signal=="(1-alpha)*10"){
          sd.sig <- (1-alpha)*10
        }else if(sd.signal=="1.5-alpha"){
          sd.sig <- 1.5-alpha
        }else if(is.numeric(sd.signal)){
          sd.sig <- sd.signal
        }else if(is.character(sd.signal)){
          error("Your term for the signal standard deviation is not yet supported.")
        }
        
        # cutoff score for alpha/2 quantile
        cut.alpha <- -qnorm(alpha/2, sd=sd.bg)
        # mean of the regulated distribution: shifted by cutoff score and beta-quantile
        mean.sig <- cut.alpha - qnorm(beta, sd=sd.sig)
        
        # Covariance matrix
        sigma.bg <- matrix(c(sd.bg^2,     sd.bg^2*rho,
                             sd.bg^2*rho, sd.bg^2     ), 2)
        sigma.both <- matrix(c(sd.sig^2,     sd.sig^2*rho,
                               sd.sig^2*rho, sd.sig^2    ), 2)
        sigma.spec1 <- matrix(c(sd.sig^2,         sd.sig*sd.bg*rho,
                                sd.sig*sd.bg*rho, sd.bg^2          ), 2)
        sigma.spec2 <- matrix(c(sd.bg^2,          sd.bg*sd.sig*rho,
                                sd.bg*sd.sig*rho, sd.sig^2         ), 2)
        
        tries <- 0
        alpha1 <- 2
        beta1 <- 2
        alpha2 <- 2
        beta2 <- 2
        stats.save <- 4.1
        
        errors <- list(errors=list(),
                       warnings=list())
        
        if(alpha<1e-4){
          diff.a <- 10*alpha
        }else if(alpha>=1e-4 & alpha<0.01){
          diff.a <- 2*alpha
        }else if(alpha>=0.01 & alpha<0.1){
          diff.a <- alpha/3
        }else if(alpha>=0.1){
          diff.a <- 0.01
        }
        if(beta<1e-4){
          diff.b <- 10*beta
        }else if(beta>=1e-4 & beta<0.01){
          diff.b <- 2*beta
        }else if(beta>=0.01 & beta<0.1){
          diff.b <- beta/3
        }else if(beta>=0.1){
          diff.b <- 0.01
        }
        
        while((abs(alpha1-alpha)>diff.a | abs(beta1-beta)>diff.b |
               abs(alpha2-alpha)>diff.a | abs(beta2-beta)>diff.b) & tries<50){
          tries <- tries+1
          
          
          ### setup data list ------------------------------------------------------
          sim.data <- list(spec1 = data.frame(id = hidden2,
                                              statistic = NA,
                                              pvalue = NA,
                                              significant = NA,
                                              direction = NA,
                                              row.names = hidden2,
                                              stringsAsFactors = F),
                           spec2 = data.frame(id = hidden2,
                                              statistic = NA,
                                              pvalue = NA,
                                              significant = NA,
                                              direction = NA,
                                              row.names = hidden2,
                                              stringsAsFactors = F),
                           meta.data = list(active.terms=active.terms,
                                            act.terms = act.terms,
                                            active.genes=active.genes,
                                            act.genes.spec1=act.genes.spec1,
                                            act.genes.spec2=act.genes.spec2,
                                            alpha = alpha,
                                            beta = beta,
                                            iteration = k,
                                            score.cutoff = cut.alpha,
                                            mean.sig = mean.sig,
                                            sd.bg = sd.bg,
                                            sd.sig = sd.sig))
          
          ### sample background ----------------------------------------------------
          # from MASS package
          bg <- data.frame(mvrnorm(length(hidden2),
                                   mu = c(0,0),
                                   Sigma = sigma.bg))
          colnames(bg) <- c("bg.spec1", "bg.spec2")
          rownames(bg) <- hidden2
          
          # initialize with background
          sim.data[["spec1"]][spec1, "statistic"] <- bg[spec1, "bg.spec1"]
          sim.data[["spec2"]][spec2, "statistic"] <- bg[spec2, "bg.spec2"]
          
          
          ### sample signal ------------------------------------------------------
          stats <- 4
          y = 0
          while((abs(alpha1-alpha)>diff.a | abs(beta1-beta)>diff.b |
                 abs(alpha2-alpha)>diff.a | abs(beta2-beta)>diff.b) & y<100){
            # if(y>0){print(paste0("Resampling necessary (", y, ")!"))}
            
            y <- y+1
            sim.data.temp <- sim.data
            
            # sample signal
            # 3 categories: active in spec1, active in spec2, active in both
            sig.both <- intersect(act.genes.spec1$gene, act.genes.spec2$gene)
            signal.both <- NULL
            
            if(length(sig.both)>0){
              if(length(sig.both)==1){
                signal.both <- data.frame(matrix(mvrnorm(length(sig.both),
                                                         mu = c(mean.sig, mean.sig),
                                                         Sigma = sigma.both),
                                                 ncol = 2))
              }else{
                signal.both <- data.frame(mvrnorm(length(sig.both),
                                                  mu = c(mean.sig, mean.sig),
                                                  Sigma = sigma.both))
              }
              colnames(signal.both) <- c("sig.spec1", "sig.spec2")
              rownames(signal.both) <- sig.both
            }
            
            sig.spec1 <- act.genes.spec1$gene[!(act.genes.spec1$gene %in% sig.both)]
            signal.spec1 <- NULL
            if(length(sig.spec1)>0){
              if(length(sig.spec1)==1){
                signal.spec1 <- data.frame(matrix(mvrnorm(length(sig.spec1),
                                                          mu = c(mean.sig, 0),
                                                          Sigma = sigma.spec1),
                                                  ncol = 2))
              }else{
                signal.spec1 <- data.frame(mvrnorm(length(sig.spec1),
                                                   mu = c(mean.sig, 0),
                                                   Sigma = sigma.spec1))
                colnames(signal.spec1) <- c("sig.spec1", "sig.spec2")
                rownames(signal.spec1) <- sig.spec1
              }
            }
            
            sig.spec2 <- act.genes.spec2$gene[!(act.genes.spec2$gene %in% sig.both)]
            signal.spec2 <- NULL
            if(length(sig.spec2)>0){
              if(length(sig.spec2)==1){
                signal.spec2 <- data.frame(matrix(mvrnorm(length(sig.spec2),
                                                          mu = c(0, mean.sig),
                                                          Sigma = sigma.spec2),
                                                  ncol = 2))
              }else{
                signal.spec2 <- data.frame(mvrnorm(length(sig.spec2),
                                                   mu = c(0, mean.sig),
                                                   Sigma = sigma.spec2))
                colnames(signal.spec2) <- c("sig.spec1", "sig.spec2")
                rownames(signal.spec2) <- sig.spec2
              }
            }
            
            # save results, in case it does not get better
            if(stats<stats.save){
              stats.save <- stats
              sim.data.save <- sim.data
              sig.both.save <- sig.both
              sig.spec1.save <- sig.spec1
              sig.spec2.save <- sig.spec2
              alpha1.save <- alpha1
              beta1.save <- beta1
              alpha2.save <- alpha2
              beta2.save <- beta2
            }
            
            # overwrite bg with signal
            if(length(sig.both)!=0){
              sim.data[["spec1"]][intersect(spec1, sig.both), "statistic"] <-
                active.genes[intersect(spec1, sig.both), "sign"]*signal.both[intersect(spec1, sig.both), "sig.spec1"]
              sim.data[["spec2"]][intersect(spec2, sig.both), "statistic"] <-
                active.genes[intersect(spec2, sig.both), "sign"]*signal.both[intersect(spec2, sig.both), "sig.spec2"]
            }
            
            sim.data[["spec1"]][intersect(spec1, sig.spec1), "statistic"] <-
              active.genes[intersect(spec1, sig.spec1), "sign"]*signal.spec1[intersect(spec1, sig.spec1), "sig.spec1"]
            sim.data[["spec2"]][intersect(spec2, sig.spec1), "statistic"] <-
              active.genes[intersect(spec2, sig.spec1), "sign"]*signal.spec1[intersect(spec2, sig.spec1), "sig.spec2"]
            
            if(length(sig.spec2)!=0){
              sim.data[["spec1"]][intersect(spec1, sig.spec2), "statistic"] <-
                active.genes[intersect(spec1, sig.spec2), "sign"]*signal.spec2[intersect(spec1, sig.spec2), "sig.spec1"]
              sim.data[["spec2"]][intersect(spec2, sig.spec2), "statistic"] <-
                active.genes[intersect(spec2, sig.spec2), "sign"]*signal.spec2[intersect(spec2, sig.spec2), "sig.spec2"]
            }
            
            # Note signs and significance
            sim.data[["spec1"]][, "sign"] <- sign(sim.data[["spec1"]][, "statistic"])
            sim.data[["spec1"]][, "significant"] <-
              (abs(sim.data[["spec1"]][, "statistic"]) > sim.data[["meta.data"]][["score.cutoff"]])
            
            sim.data[["spec2"]][, "sign"] <- sign(sim.data[["spec2"]][, "statistic"])
            sim.data[["spec2"]][, "significant"] <-
              (abs(sim.data[["spec2"]][, "statistic"]) > sim.data[["meta.data"]][["score.cutoff"]])
            
            contig1 <- table(sim.data[["spec1"]]$significant>0,
                             rownames(sim.data[["spec1"]]) %in% c(sig.spec1, sig.both))
            
            if(dim(contig1)[1]==2){ # check for full table
              alpha1 <- contig1[2, 1]/(contig1[2, 1] + contig1[1, 1])
              beta1 <- contig1[1, 2]/(contig1[2, 2] + contig1[1, 2])
            }else{ # attach 0s
              alpha1 <- as.numeric(as.logical(rownames(contig1)))
              beta1 <- 1-as.numeric(as.logical(rownames(contig1)))
            }
            
            contig2 <- table(sim.data[["spec2"]]$significant>0,
                             rownames(sim.data[["spec2"]]) %in% c(sig.spec2, sig.both))
            if(dim(contig2)[1]==2){ # check for full table
              alpha2 <- contig2[2, 1]/(contig2[2, 1] + contig2[1, 1])
              beta2 <- contig2[1, 2]/(contig2[2, 2] + contig2[1, 2])
            }else{ # attach 0s
              alpha2 <- as.numeric(as.logical(rownames(contig2)))
              beta2 <- 1-as.numeric(as.logical(rownames(contig2)))
            }
            
            stats.new <- sum(abs(alpha1-alpha),
                             abs(beta1-beta),
                             abs(alpha2-alpha),
                             abs(beta2-beta))
            if(stats.new>stats &
               (abs(alpha1-alpha)>diff.a | abs(beta1-beta)>diff.b |
                abs(alpha2-alpha)>diff.a | abs(beta2-beta)>diff.b)){
              sim.data <- sim.data.temp
            }else{
              stats <- stats.new
            }
          }
          
          print(paste0("Total difference to desired error rates: ",
                       stats,
                       " after resampling ", y, " times."))
          
          if(stats>stats.save &
             (abs(alpha1-alpha)>diff.a | abs(beta1-beta)>diff.b |
              abs(alpha2-alpha)>diff.a | abs(beta2-beta)>diff.b)){
            print("Falling back to saved results")
            # save results, in case it does not get better
            stats <- stats.save
            sim.data <- sim.data.save
            sig.both <- sig.both.save
            sig.spec1 <- sig.spec1.save
            sig.spec2 <- sig.spec2.save
            alpha1 <- alpha1.save
            beta1 <- beta1.save
            alpha2 <- alpha2.save
            beta2 <- beta2.save
          }else if(stats<stats.save &
                   (abs(alpha1-alpha)>diff.a | abs(beta1-beta)>diff.b |
                    abs(alpha2-alpha)>diff.a | abs(beta2-beta)>diff.b)){
            stats.save <- stats
            sim.data.save <- sim.data
            sig.both.save <- sig.both
            sig.spec1.save <- sig.spec1
            sig.spec2.save <- sig.spec2
            alpha1.save <- alpha1
            beta1.save <- beta1
            alpha2.save <- alpha2
            beta2.save <- beta2
            
            print("Update fallback!")
          }
        }
        
        if((abs(alpha1-alpha)>diff.a | abs(beta1-beta)>diff.b |
            abs(alpha2-alpha)>diff.a | abs(beta2-beta)>diff.b)){
          warning(paste0("Sampling unsuccessful: iteration ", k,
                         ", alpha ", alpha, ", beta ", beta))
          print(paste0("FPR = ", alpha1,  "(", alpha, ")",
                       ", FNR = ", beta1,  "(", beta, ")"))
          print(paste0("Total difference to desired error rates: ",
                       stats,
                       " after resampling ", y, " times."))
        }
        
        ## Diagnostics: ----------------------------------------------------------
        
        #.........................................................................
        # #                             Truth: PW regulated
        # #                              FALSE    |   TRUE
        # #                        _______________|_______________
        # #               FALSE   |     TN            FN
        # # prediction   _________|
        # # PW regulated          |
        # #               TRUE    |     FP            TP
        # #
        #
        #
        # pnorm(qnorm(alpha/2, sd=sd.bg), mean = -(qnorm(alpha/2, sd=sd.bg) + qnorm(beta, sd=sd.sig)), sd=sd.sig)
        
        # # False positive rate FPR (alpha) = type I error = FP / (FP + TN)
        #                                                  = FP / (all negative)
        # alpha
        # contig[2, 1]/(contig[2, 1] + contig[1, 1])
        #
        # # False negative rate FNR (beta) = type II error = FN / (TP + FN)
        #                                                  = FN / (all positive)
        # beta
        # contig[1, 2]/(contig[2, 2] + contig[1, 2])
        #
        
        ### Species 1 ------------------------------------------------------------
        contig <- table(sim.data[["spec1"]]$significant>0,
                        rownames(sim.data[["spec1"]]) %in% c(sig.spec1, sig.both))
        contig
        
        if(dim(contig)[1]==2){ # check for full table
          alpha
          alpha.sim <- contig[2, 1]/(contig[2, 1] + contig[1, 1])
          
          beta
          beta.sim <- contig[1, 2]/(contig[2, 2] + contig[1, 2])
          
          qc.spec1 <- rbind(qc.spec1,
                            c(k,
                              alpha,
                              alpha.sim,
                              beta,
                              beta.sim,
                              abs(alpha.sim-alpha)>diff.a,
                              abs(beta.sim-beta)>diff.b))
        }else{ # attach 0s
          qc.spec1 <- rbind(qc.spec1,
                            c(k,
                              alpha,
                              as.numeric(as.logical(rownames(contig))),
                              beta,
                              1-as.numeric(as.logical(rownames(contig))),
                              abs(as.numeric(as.logical(rownames(contig)))-alpha)>diff.a,
                              abs(1-as.numeric(as.logical(rownames(contig)))-beta)>diff.b))
          alpha
          beta
        }
        
        ### Species 2 ------------------------------------------------------------
        contig <- table(sim.data[["spec2"]]$significant>0,
                        rownames(sim.data[["spec2"]]) %in% c(sig.spec2, sig.both))
        contig
        
        if(dim(contig)[1]==2){ # check for full table
          alpha
          alpha.sim <- contig[2, 1]/(contig[2, 1] + contig[1, 1])
          
          beta
          beta.sim <- contig[1, 2]/(contig[2, 2] + contig[1, 2])
          
          qc.spec2 <- rbind(qc.spec2,
                            c(k,
                              alpha,
                              alpha.sim,
                              beta,
                              beta.sim,
                              abs(alpha.sim-alpha)>diff.a,
                              abs(beta.sim-beta)>diff.b))
        }else{ # attach 0s
          qc.spec2 <- rbind(qc.spec2,
                            c(k,
                              alpha,
                              as.numeric(as.logical(rownames(contig))),
                              beta,
                              1-as.numeric(as.logical(rownames(contig))),
                              abs(as.numeric(as.logical(rownames(contig)))-alpha)>diff.a,
                              abs(1-as.numeric(as.logical(rownames(contig)))-beta)>diff.b))
          alpha
          beta
        }
        
        print(paste0("Iteration ", k,
                     ", combination alpha = ", alpha,
                     ", beta = ", beta,
                     " done ( ", l, " / ", dim(error.rates)[1], " )"))
        
        saveRDS(sim.data,
                file=paste0(outdir, "/", label, "/sum_stats/",
                            "sim_data_", k, "_", alpha, "_", beta, ".RDS"))
        
      }
      
      sim_all <- list()
      
      for (l in 1:dim(error.rates)[1]){
        alpha <- error.rates[l, "alpha"]
        beta <- error.rates[l, "beta"]
        sim.data <- readRDS(file=paste0(outdir, "/", label, "/sum_stats/",
                                        "sim_data_", k, "_", alpha, "_", beta, ".RDS"))
        sim_all[[paste0(alpha, "_", beta)]] <- sim.data
      }
      saveRDS(sim_all,
              file=paste0(outdir, "/", label, "/sum_stats/",
                          "sim_data_all_", k, ".RDS"))
      print(paste0("Data collected for iteration ", k,
                   " done ( ", k, " / ", length(iterations), " )"))
      
      for (l in 1:dim(error.rates)[1]){
        alpha <- error.rates[l, "alpha"]
        beta <- error.rates[l, "beta"]
        file.remove(paste0(outdir, "/", label, "/sum_stats/",
                           "sim_data_", k, "_", alpha, "_", beta, ".RDS"))
      }
      print(paste0("Unneeded data deleted for iteration ", k,
                   " ( ", k, " / ", length(iterations), " )"))
      
      print(paste0("Iteration ", k, " done ( ", k, " / ", length(iterations), " ) for scenario ", label))
      print("Active pathways:")
      print.data.frame(sim.data$meta.data$active.terms,
                       row.names = F)
    }
    saveRDS(qc.spec1,
            file=paste0(outdir, "/", label, "/QC_stats/",
                        "stats_species1.RDS"))
    saveRDS(qc.spec2,
            file=paste0(outdir, "/", label, "/QC_stats/",
                        "stats_species2.RDS"))
    
    qc.spec1 <- readRDS(paste0(outdir, "/", label,
                               "/QC_stats/stats_species1.RDS"))
    colnames(qc.spec1) <- c("iteration", "alpha", "alpha_sim",
                            "beta", "beta_sim",
                            "fail_alpha", "fail_beta")
    
    qc.spec2 <- readRDS(paste0(outdir, "/", label,
                               "/QC_stats/stats_species2.RDS"))
    colnames(qc.spec2) <- c("iteration", "alpha", "alpha_sim",
                            "beta", "beta_sim",
                            "fail_alpha", "fail_beta")
    
    qc.spec1$dif.alpha <- qc.spec1$alpha-qc.spec1$alpha_sim
    qc.spec1$dif.beta <- qc.spec1$beta-qc.spec1$beta_sim
    qc.spec2$dif.alpha <- qc.spec2$alpha-qc.spec2$alpha_sim
    qc.spec2$dif.beta <- qc.spec2$beta-qc.spec2$beta_sim
    
    qc.spec1 <- qc.spec1[order(abs(qc.spec1$dif.beta), decreasing = T), ]
    qc.spec2 <- qc.spec2[order(abs(qc.spec2$dif.beta), decreasing = T), ]
    
    saveRDS(list(qc.spec1=qc.spec1, qc.spec2=qc.spec2),
            file=paste0(outdir, "/", label,
                        "/QC_stats/stats_all.RDS"))
    
    range <- max(c(abs(qc.spec1$dif.alpha), abs(qc.spec1$dif.beta),
                   abs(qc.spec2$dif.alpha), abs(qc.spec2$dif.beta)))
    print(paste0("Maximal deviation from desired error rate: ", range))
    parameters$max.error <- range
    
    pdf(paste0(outdir, "/", label,
               "/QC_stats/stats_all.pdf"),
        height=15, width=20)
    par(mfrow=c(2,2))
    hist(qc.spec1$dif.alpha, breaks = 20,
         main = "FPR alpha, Species 1",
         xlim = c(-range, range))
    hist(qc.spec1$dif.beta, breaks = 20,
         main = "FNR beta, Species 1",
         xlim = c(-range, range))
    hist(qc.spec2$dif.alpha, breaks = 20,
         main = "FPR alpha, Species 2",
         xlim = c(-range, range))
    hist(qc.spec2$dif.beta, breaks = 20,
         main = "FNR beta, Species 2",
         xlim = c(-range, range))
    par(mfrow=c(1,1))
    dev.off()
    
    qc1 <- qc.spec1[, c("alpha", "beta", "dif.alpha", "dif.beta")]
    qc1$dif.alpha <- abs(qc1$dif.alpha)
    qc1$dif.beta <- abs(qc1$dif.beta)
    qc1 <- aggregate(. ~alpha+beta, data=qc1, sum, na.rm=TRUE)
    qc2 <- qc.spec2[, c("alpha", "beta", "dif.alpha", "dif.beta")]
    qc2$dif.alpha <- abs(qc2$dif.alpha)
    qc2$dif.beta <- abs(qc2$dif.beta)
    qc2 <- aggregate(. ~alpha+beta, data=qc2, sum, na.rm=TRUE)
    
    library(ggplot2)
    library(ggpubr)
    pdf(paste0(outdir, "/", label,
               "/QC_stats/stats_all_error.pdf"),
        height=15, width=20)
    print(ggarrange(
      ggplot(qc1, aes(x=alpha, y=beta,
                      fill=dif.alpha)) +
        geom_tile() +
        ggtitle("Species 1, alpha"),
      ggplot(qc1, aes(x=alpha, y=beta,
                      fill=dif.beta)) +
        geom_tile() +
        ggtitle("Species 1, beta"),
      ggplot(qc1, aes(x=alpha, y=beta,
                      fill=dif.alpha+dif.beta)) +
        geom_tile() +
        ggtitle("Species 1, alpha+beta"),
      ggplot(qc2, aes(x=alpha, y=beta,
                      fill=dif.alpha)) +
        geom_tile() +
        ggtitle("Species 2, alpha"),
      ggplot(qc2, aes(x=alpha, y=beta,
                      fill=dif.beta)) +
        geom_tile() +
        ggtitle("Species 2, beta"),
      ggplot(qc2, aes(x=alpha, y=beta,
                      fill=dif.alpha+dif.beta)) +
        geom_tile() +
        ggtitle("Species 2, alpha+beta"),
      nrow = 2, ncol = 3,
      common.legend = F
    ))
    dev.off()
    
    print("Simulation complete!")
    
  })
  
  print(time)
  parameters$time.of.simulation <- paste0("User: ", time[1],
                                          ", System: ", time[2],
                                          ", Elapsed: ", time[3])
  
  # write parameters to file
  write.table(data.frame(t(parameters)),
              file=paste0(outdir, "/", label, "/parameter_info.txt"),
              row.names = T, col.names = F, quote = F, sep = "\t")
  
  return(errors)
  
}
