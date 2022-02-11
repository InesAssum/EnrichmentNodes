#' Extension to the waitForJobs function of the BatchJobs package which shows
#' some strange behaviour when waiting for jobs (database locked)
#' so we need to make it extra failsafe.
#'
#' @family eqtl functions
#' @title basic eQTL interface
#' @param reg batchtools registry
#' @param waittime time to wait before updating job status
#' @param nretry number of time to retry getting the job status before throwing
#'        an error
#'
#' @imports batchtools
#' 
#' 

# rfun = do.stuff
# idx = list(1:3)
# more.args = list()
# name = "test"
# dir = "tmp_dir_test"
# clean.up = F
# resources = list(time="00:05:00", memory="200M")

myWaitForJobs <- function(ids, reg, waittime=3, nretry=100) {
  require(batchtools)
  success = FALSE
  while (nretry > 0 && !success) {
    status = tryCatch({
      while (TRUE) {
        status = getStatus(ids = ids, reg = reg)
        print(status)
        if (status$done + status$error + status$expired == status$defined) {
          cat("done\n")
          return(list(success=TRUE, nretry=nretry))
        }
        Sys.sleep(waittime)
      }
      return(list(success=FALSE, nretry=nretry))
    }, error=function(e) {
      cat("Error while waiting for jobs:\n")
      print(e)
      cat("\nnumber of retries left: ", nretry - 1, "\n")
      Sys.sleep(waittime + runif(1, 0, 3))
      return(list(success=FALSE, nretry=nretry - 1))
    })
    success = status$success
    nretry = status$nretry
    cat("success after the tryCatch block:", success, "\n")
    cat("nretry after the tryCatch block:", nretry, "\n")
  }
  
  
  if (!success) {
    err.msg = paste("Error during batch processing in registry")
    save(envir=sys.frame(), list=ls(envir=sys.frame()), file=file.path(dir, "error_image.RData"))
    stop(err.msg)
  }
}


run.batchtools <- function(rfun, idx, more.args, name=NULL, dir, clean.up=TRUE,
                           resources=list(), n.chunks=NULL, reuse.registry=FALSE) {
  require(batchtools)
  
  if (reuse.registry) {
    reg <- loadRegistry(file.dir = dir, writeable = TRUE)
    reg <- getDefaultRegistry()
  } else {
    reg = makeRegistry(file.dir=dir)
  }
  print("Registry created")
  ids <- batchMap(fun = rfun,
                  idx,
           more.args=more.args,
           reg = reg)
  if (length(name)==1) {
    setJobNames(ids, paste0(name, idx), reg = reg)
  } else if (length(name)==length(idx)) {
    setJobNames(ids, name, reg = reg)
  } else {
    setJobNames(ids, paste0(gsub("tmp_dir_", "", dir), idx), reg = reg)
  }
 
  
  ## jobs can be packed together in chunks if there are too many
  library(data.table)
  if (!is.null(n.chunks)) {
    ids[, chunk := chunk(job.id, n.chunks = n.chunks)]
  }
  print("checked for chunks")
  
  ## job delay to prevent concurrent access to the database by too many jobs
  print("Let's submit jobs")
  submitJobs(ids = ids, reg = reg, resources=resources)
  print("jobs submitted")
  Sys.sleep(20)
  
  print("Let's wait!")
  
  ## also a custom wait function that is more error tolerant with many jobs
  myWaitForJobs(ids=ids, reg=reg, waittime=30, nretry=100)
  Sys.sleep(20)
  res = reduceResultsList(reg=reg)
  job = getJobTable(reg=reg)
  print(head(getJobTable(reg=reg), 10))
  if (clean.up) {
    removeRegistry(wait=0, reg=reg)
  }
  return(list(res=res, job=job))
}


## Mini example:
#
# library(batchtools)
# rfun <- function(i, filename) {
#   # do some very important stuff
#   print(filename)
#   Sys.sleep(1200)
#   print(i)
#   return(i)
# }

# test_batch <- run.batchtools(rfun, idx = 1:15, more.args = list(filename="123"),
#                              name = "sim", dir = "tmp_dir_simulation3",
#                              clean.up=F,
#                              resources=list(partition="icb_cpu",
#                                             memory="200M",
#                                             ncpus = 1,
#                                             measure.memory = TRUE,
#                                             walltime="00:15:00"))
#
# test_batch <- run.batchtools(rfun, 1:10, more.args=list(filename="123"),
#                              "sim", "tmp_dir_simulation2",
#                              clean.up=F, n.chunks = 5,
#                              resources=list(partition="icb_cpu",
#                                             memory="200M",
#                                             ncpus = 1,
#                                             measure.memory = TRUE,
#                                             walltime="00:05:00"))
# 
# str(test_batch)

# load("/home/icb/ines.assum/rstudio/tmp_dir_simulation/registry.RData")
# submitJobs(reg, findErrors(reg))
# submitJobs(reg, findExpired(reg))