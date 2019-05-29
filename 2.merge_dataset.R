##################################################
## Description: Merging GBIF, Hulten and PAF datasets
## 
## Date: 2018-08-12 11:49:22
## Author: Oskar Hagen (oskar@hagen.bio)
##################################################

#set working directory
setwd("SET YOUR WORKING DIRECTORY HERE")

#source relevant functions
source("./support_functions/resolve_taxonomic_names.R")
source("./support_functions/merge_occurance_data.R")
source("./support_functions/compress_housekeeping.R")

OutputFolder <- "SELECT DESIRED OUTPUT FOLDER"

#resolve species on target datasets
sdf <- resolve.sps(Input_folder_s=c("../1_gbif/point","../2_hulten/point", "../3_paf/point"))

#store sdf for future reference
write.csv2(sdf, file=file.path(OutputFolder, "script", "sdf.csv"))

#run merging function
merge.sdf(sdf, output_folder=OutputFolder)