##################################################
## Description: Support function declaration
## This function merges species occurance data in a same dataset 
## with the same finalname based on sdf matrix. This should be used after resolving 
## automaticaly and manualy species taxonomic names
##
## Project: Distribution of Cold Adapted Plants data as points 
##
## Date: 2018-03-19 21:59:00
## Author: Oskar Hagen (oskar@hagen.bio)
##################################################


merge.sdf <- function(sdf, output_folder){
  
  #Function description
  # This function merges species with the same finalname based on sdf matrix.
  
  #Parameters description
  # sdf: contain the matrix or data.frame generated with resolve.sps (sourced with resolve_taxonomic_names.R)
  # output_foler: is the place were merged files will be created
  
  
  #merge species on output_folder based on data.frame with "path", "submittedname" and "finalname.
  sdf <- as.data.frame(sdf)
  #merge distributions based on a reference table, i.e. finalname. We need to have the location of
  #the folders
  if (file.exists(file.path(output_folder, "merged_point"))){
    stop("merged_point folder already exist! Delete manually and run data.sets function again if this is what you wannt")
  } else {
    dir.create(file.path(output_folder, "merged_point"))
  }
  
  un <- unique(sdf$finalname) # get unique species
  
  
  pb <- txtProgressBar(min = 0, max = length(un), style = 3)
  for (uni in un){ #loop over unique species
    # print(paste("working with", uni, "index=", match(uni, un)))
    presences <- NULL
    presences <- which(sdf$finalname==uni)
    mpts <- NULL
    for (pi in presences){ #loop over datasets containing the points...
      p <- read.table(file.path(sdf[pi,"path"], paste0(sdf[pi,"submittedname"],".txt")), header = T)[,1:2]
      colnames(p) <- c("X", "Y")
      # print(dim(p))
      mpts <- rbind(mpts, p)
      # remove duplicates
      mpts <- unique(mpts)
    }
    # print(dim(mpts))
    # save points...
    write.table(mpts, file = file.path(output_folder, "merged_point", paste0(uni, ".txt") ), row.names = FALSE)
    setTxtProgressBar(pb, match(uni, un))
  }# end loop over unique species....
  close(pb)
}