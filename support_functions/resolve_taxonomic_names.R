##################################################
## Description: Support function declaration
## Resolve taxonomic name and look for synonyms while deleting clear mistakes
## Project: Distribution of Cold Adapted Plants data as points
##
## Details: using Taxonomic Name Resolution Service via the Taxonstand The Plant List (TPL) 
##
## Date: 2018-03-19 21:59:00
## Author: Oskar Hagen (oskar@hagen.bio)
##################################################

library("Taxonstand")
library("utils")


resolve.sps <- function(Input_folder_s){

  #Function description
  # !this function requires internet connection!
  # It creates a matrix with the input folders and submittedname and finalname
  # a csv is created at the outputfolder
  
  #example inputs:
  # Input_folder_s <- c("../1_gbif/point", "../2_hulton/point")
  # output_folder <- "../3_gbif_hulton"
  
  # get species names
  sdf <- list()
  for (i in 1:length(Input_folder_s)){
    #check if point folder exists
    if (!file.exists(Input_folder_s[i])){
      stop("point folder does not exist, please check this dataset")
    }
    sdf[[i]] <- cbind(Input_folder_s[i], tools::file_path_sans_ext(list.files(Input_folder_s[i])))
  }
  sdf <- do.call(rbind, sdf)
  colnames(sdf) <- c("path", "submittedname")
  
  
  
  spnames <- sdf[,"submittedname"]
  # resolve
  # m <- tnrs(query = spnames)#[ , -c(5:7)]
  mtpl <- TPL(spnames, diffchar = 4, drop.lower.level=T)
  
  
  sdf <- cbind(sdf, finalname=paste(mtpl$"New.Genus", mtpl$"New.Species"))
  
  #clean sdf ------------------
  # remove things that are not only 2 elements
  elements <-   unlist(lapply(strsplit(sdf[, "finalname"], "\\s+"), length))
  sdf <- sdf[!elements!=2, ,drop=F]
  nchars <- nchar(sdf[, "finalname"], type = "chars", allowNA = FALSE, keepNA = NA)
  sdf <- sdf[nchars>2, ,drop=F]
  # remove bad characters
  # "..txt"
  sdf <- sdf[-grep(".txt", sdf[, "finalname"], fixed=T), ,drop=F]
  # sdf <- sdf[-grep("..", sdf[, "finalname"], fixed=T), ,drop=F]
  # remove NA's
  sdf <- na.omit(sdf)
  removethis <- c(" NA\\.asc", "\\d", "\\b\\. ", "^NA ", "^[a-zA-Z]{2} ", "^NA$", " NA$", "'", "\"", "/")
  
  mistakes <- grep(paste(removethis,collapse="|"), sdf[,"finalname"], value=T)
  
  print(paste0("we are removing: >", paste0(paste0(mistakes, collapse = "<  >"),"<") ))
  
  sdf <- sdf[!sdf[,"finalname"]%in%mistakes,]
  
  
  #end clean sdf---------------
  
  return(sdf)
}

#optional
#save csv
#write.csv(sdf, file.path("...../sdt.csv"), row.names=F)
