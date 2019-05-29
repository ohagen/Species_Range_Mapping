##################################################
## Description: Transform points to range raster
## simplified version for non-cluster computing
## Uses: ecoregions downloaded from 
## The Nature Conservancy (25th April 2018). Terrestrial Ecoregions of the World. https:// http://maps.tnc.org
## and climate layer downloaded from 
## GBIF.org (26th August 2018) GBIF Occurrence Download https://doi.org/10.15468/dl.1xd4qs
##
## Date: 2018-09-03 13:57:16
## Author: Oskar Hagen (oskar@hagen.bio)
##################################################

#set working directory
setwd("SET YOUR WORKING DIRECTORY HERE")

#load libraries
library(raster)
library(rgdal)
library(maptools)
library(tools)

#source relevant functions
source("./support_functions/range_from_points.R")

#load variables and define parameters
sps <- list.files('../merged_point/')
proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#load ecoregions layer
ecoregions <- readShapePoly("DOWNLOADED ECOREGIONS LAYER", proj4string = proj)
temperature <- raster('DOWNLOADED CLIMATE LAYER') 
fromhere <- 1
tohere <- length(sps)
tempfiles <- "C:/temp"
input_dir <- "../merged_point"
output_dir <- "../range_raster"

# set temp DIR...
if (!file.exists(paste0('./', tempfiles))){
  dir.create(paste0('./', tempfiles), recursive=T)
}

# loop over species
for (indexi in fromhere:tohere){
  spi <- sps[indexi]
  spi_name <- file_path_sans_ext(spi)
  
  print(paste0("working with        ", spi_name, "         [index:", indexi, "]"))
  
  #run range mapping function
  t <- range.fun(species_name = spi_name, occ_coord =  na.omit(read.table(file.path(input_dir, spi), header=T)), proj=proj, Climatic_layer=temperature, Bioreg = ecoregions, final_resolution=1, degrees_outlier=7, 
             clustered_points_outlier=1, buffer_width_point=0.5, buffer_increment_point_line=0.5, buffer_width_polygon=0.1, cut_off=0.01, method="PCM", cover_threshold=0.1,
    write_plot=F, dest_plot='./Output/Maps/', write_raster=T, dest_raster=file.path(output_dir), return_raster=T, overwrite=T, dir_temp=tempfiles, dir_log='./log/')
  
  print(paste0("COMPLETED: [ ",round(indexi/length(sps),2)*100, "% ]"))
}

#delete temp files
unlink(paste0('./', tempfiles), recursive=T)
