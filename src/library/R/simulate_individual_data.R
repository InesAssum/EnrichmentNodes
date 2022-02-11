#' #############################################################################
#  Simulate individual level transcript and protein expression  ################
#' #############################################################################
#'
#' Simulate individual level transcript and protein expression
#'
#' @param sim.data simulated data containing summary statistics and gene coverage
#' @param dir path to out directory (string) (where to save files)
#'
#' @param rho correlation between omics
#' @param N.iter number of replicates/scenarios
#' @param N.ind number of individuals to simulate
#' @param frac fraction of individuals that are cases
#'
#' @author Ines Assum

simulate.ind.data <- function(dir, sim.data,
                              rho=0, N.iter=1:100, N.ind=100, frac=0.3,
                              seed=NULL, mu=NULL, sigma=NULL){
  require(PathwayGuidedRF)
  require(Umpire)
  
  print(paste0("Computations starting at ", Sys.time()))
  
  if(!file.exists(paste0(dir, "/sim_data_mvn.RDS"))){
    if(is.null(seed)){
      seed <- as.integer(runif(1)*2e9)
    }
    set.seed(seed)
    print(paste0("Current seed is: ", seed))
    
    #sim.data <- readRDS(paste0(dir, "/sim_data_1.RDS"))[[1]][["1e-04_1e-04"]]
    genes <- c(paste0(names(sim.data[["mRNA"]]), "_trans"),
               paste0(names(sim.data[["prot"]]), "_prot"))
    sigma <- matrix(0, nrow = length(genes), ncol = length(genes),
                    dimnames = list(genes, genes))
    if(rho>0){
      ov <- intersect(names(sim.data[["mRNA"]]), names(sim.data[["prot"]]))
      for(x in ov){
        sigma[paste0(x, "_prot"),paste0(x, "_trans")] <- 1
        sigma[paste0(x, "_trans"),paste0(x, "_prot")] <- 1
      }
      sigma <- rho*sigma
    }
    diag(sigma) <- 1
    
    data <- list()
    if(is.null(mu)){mu <- rep(0, length(genes))}
    data[["mvn"]] <- MVN(mu, sigma)
    
    no.ctrl <- round(N.ind*(1-frac))
    no.case <- round(N.ind*frac)
    data[["y"]] = c(rep(0, no.ctrl), rep(1, no.case))
    
    data[["rand"]] <- rand(data[["mvn"]], N.ind)
    rownames(data[["rand"]]) <- genes
    
    data[["SessionInfo"]] <- sessionInfo()
    
    saveRDS(data, file=paste0(dir, "/sim_data_mvn.RDS"))
    #saveRDS(mvn, file="mvn_temp.RDS")
  } else {
    data <- readRDS(paste0(dir, "/sim_data_mvn.RDS"))
    seed <- data[["seed"]]
    set.seed(seed)
    print(paste0("Current seed is: ", seed))
  }
    
  fpr <- (0:10)/10
  fpr[c(1,11)] <- c(1e-4, 1-1e-4)
  fnr <- (0:10)/10
  fnr[c(1,11)] <- c(1e-4, 1-1e-4)
  
  for (k in N.iter){
    
    start1 <- Sys.time()
    
    sim.data <- readRDS(paste0(dir, "/sim_data_", k, ".RDS"))[[paste(k)]]
    dir.create(paste0(dir, "/ind_data/", k), recursive = T, showWarnings = F)
    
    l <- 0
    istart <- 1
    
    for (i in istart: length(fpr)) {
      if (istart > 1) {istart <- 1}
      print(paste0("Iteration ", k-min(N.iter)+1, " of ", max(N.iter)-min(N.iter)+1, " running: ", l, " / ", length(sim.data)))
      for (j in 1: length(fnr)) {
        # i=1
        # j=1
    
        alpha <- fpr[i]
        beta <- fnr[j]
        
        outfile <- paste0(dir, "/ind_data/", k, "/sim_data_ind_", alpha, "_", beta, ".RDS")
        
        if(!file.exists(outfile) & !is.null(sim.data[[paste0(alpha, "_", beta)]][["mRNA"]])){
          
          l <- l+1
          
          rank.mRNA <- sim.data[[paste0(alpha, "_", beta)]][["mRNA"]]
          rank.prot <- sim.data[[paste0(alpha, "_", beta)]][["prot"]]
          
          sd.bg <- 1+100*beta
          sd.sig <- (1-alpha)*10
          cut.alpha <- -qnorm(alpha/2, sd=sd.bg)
          
          full.data <- list()
          
          full.data[["mRNA"]] <- as.matrix(data[["rand"]][paste0(names(rank.mRNA), "_trans"), ])
          full.data[["mRNA"]][, data[["y"]]==1] <- full.data[["mRNA"]][, data[["y"]]==1] + 5*as.numeric(scale(rank.mRNA, center = F))/cut.alpha
          rownames(full.data[["mRNA"]]) <- names(rank.mRNA)
          full.data[["mRNA"]] <- data.frame(y=data[["y"]],
                                            t(full.data[["mRNA"]]))
          
          full.data[["prot"]] <- as.matrix(data[["rand"]][paste0(names(rank.prot), "_prot"), ])
          full.data[["prot"]][, data[["y"]]==1] <- full.data[["prot"]][, data[["y"]]==1] + 5*as.numeric(scale(rank.prot,
                                                                                          center = F))/cut.alpha
          rownames(full.data[["prot"]]) <- names(rank.prot)
          full.data[["prot"]] <- data.frame(y=data[["y"]],
                                            t(full.data[["prot"]]))
          
          full.data[["seed"]] <- seed
          
          saveRDS(full.data,
                  file = outfile)
          
        }
      }
    }
    end1 <- Sys.time()
    print(paste0("Iteration ", k, " took ", end1-start1))
  }
}



