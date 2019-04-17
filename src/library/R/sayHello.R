#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (!(length(args)==2)) {
  stop("Two argument must be supplied ((input file).n)string)", call.=FALSE)
}

sayHello <- function(){
  writeLines(args[1], args[2])
}
sayHello()
