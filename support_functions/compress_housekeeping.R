##################################################
## Description: Support function declaration
## Compress point data and delete original folder
##
## Project: Distribution of Cold Adapted Plants data as points 
##
## Date: 2018-03-19 21:59:00
## Author: Oskar Hagen (oskar@hagen.bio)
##################################################


compress.point <- function(Input_folder){
  # compress point and delete original folder!!! attention!!
  # input folder called "point" is assumed
  bkpdir <- getwd()
  setwd(file.path(Input_folder, "point")) 
  files2zip <- dir('.', full.names = TRUE)
  zip(zipfile = "../point", files = files2zip)
  # warnning : DELETE ALL files inside folder POINT
  unlink("../point/*", recursive=T)
  file.rename(from = "../point.zip",  to = "./point.zip")
  #set back wd
  setwd(bkpdir)
  return(0)
}