#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#' HelloWorld script
#' 
#' @param args[1] String/message to print.
#' @param args[2] Filename to save results.
#' 
#' @author Ines Assum

sayHello <- function(){
  writeLines(args[1], args[2])
  print(args[1])
}
sayHello()





