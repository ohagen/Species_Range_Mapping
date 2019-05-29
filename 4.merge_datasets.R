##################################################
## Description: Merging PanArtic with (GBIF already with Hulten)
## this script do not require us to have PanArtic already taxonomicaly solved...
##
## Date: 2018-01-16 16:18:22
## Author: Oskar Hagen (oskar@hagen.bio)
##################################################

#load libraries
library("raster")
library("rgdal")
library("maptools")

#set R variables
rasterOptions(tmpdir="S:/temp")

#define top directory
topdir <- "S:/ohagen/ColdPlants/scripts/RetriveDistData/5_gbif_hulten_panarticflora"

h_dir <- file.path(topdir,"../3_gbif_hulten")
p_dir <- file.path(topdir, "../4_panarticflora")
output_dir <- file.path(topdir, "ranges_raster")

#defina variables
proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
resolution <- "1d" #1 degree final resolution for distribution rasters

# prepare masks
# "LOAD HERE TEMPERATURE LAYER"
#DO NOT RUN
mask <- "mask containing 1 on land and NA on ocean"



#### start internal functions ####
savetodisk <- function(raster, sp_name, output_dir=output_dir){
  #saves to disk without considering land areas
  raster[][raster[]==0] <- NA # remove land 
  #crop final raster extent...
  raster <- crop(raster, extent(-180, 180, 10, 90))
  writeRaster(raster, filename=paste( output_dir, "/", sp_name,".asc",sep=""),overwrite=TRUE)
  print(paste("SAVED TO DISK", sp_name, "resolution:", res(raster)[1]))
}

savetodiskmaskrasta <- function(raster, mask, sp_name, output_dir="../mask_raster"){
  #saves to disk considering a mask
  raster <- raster*mask
  writeRaster(raster, filename=paste( output_dir, "/", sp_name,".asc",sep=""),overwrite=TRUE)
  print(paste("SAVED TO DISK", sp_name, "resolution:", res(raster)[1]))
}
#### end internal functions ####


#read dummy rasta
dummyrasta <- raster(file.path(topdir, "dummyrasta.asc")) #here you can load a raster with same extension without information

#read PANArticFlora
pandata <- readRDS(paste0(p_dir, "/scripts/PANARTICDATAHERE"))

#load resolved PANArticFlora
pandata <- readRDS(file.path(topdir, "FILE GENERATED WITH resolve_taxonomic_names.R and checked/modified manually"))

#deal with specificities from PANArcticFlora
zones <- readShapePoly(paste0(p_dir, "/scripts/zones/PanArticZones"))
projection(zones) <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"
zones <-spTransform(zones, proj)

#create output directory
output_dir_res <- file.path(output_dir, resolution)
if (!dir.exists(output_dir_res)){
  dir.create(output_dir_res)
}

# 1 read GBIFhulten data
sh <- list.files(file.path(h_dir,"range_raster" ) )

#createlist
sdf <- cbind("hulten", gsub(".asc", "", sh, fix=T))
sdf <- rbind(sdf, cbind("pan", unique(pandata$X.1)))
colnames(sdf) <- c("path", "name")
sdf <- as.data.frame(sdf, stringsAsFactors=F)

# read already finished runs
alreadyrun <- list.files(file.path(output_dir, resolution ) )
alreadyrun <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(alreadyrun))

un <- unique(sdf$name)

#loop over unique species
for (lhi in 1:length(un)){
  
  spi <- un[lhi]
  
  print(paste0("lhi: ", lhi, "  | species:  [ ", spi, " ]   [", round((lhi/length(un))*100, 0), "%] - ", format(Sys.time(), '%H:%M:%S') ))
  
  rasta <- dummyrasta
  
  prob_loading <- F
  
  sc <- sdf[sdf$name==spi,]
  

  if ("hulten"%in%sc$path){ #if in hulten
     
     rasta <- tryCatch(raster(file.path(h_dir, "range_raster", paste0(spi, ".asc") )), 
                       error=function(e){
                                         print("prob. loading rasta")
                                         prob_loading<<-T
                                         write.table("error loading...", file=file.path(topdir, "prob_loading", paste0(spi,".txt")), col.names=F, row.names=F)
                                         })
     if (prob_loading){
       rasta <- dummyrasta
     }
  }
  
  #quality check plot
  ## jpeg(filename=file.path(topdir,"plot", paste0(spi, ".jpeg") ), width = 580, height = 580)
  ## 
  ## par(mfrow=c(2,1))
  ## plot(crop(rasta, extent(-180, 180, 10, 90)), legend=FALSE, axes=FALSE, box=FALSE)
  ## legend(x=-150, y=95, horiz=T, legend=c("Hulten", "Panarticflora"),bty='n', pch=15, col=c("darkgreen", "blue"))
  
  if ("pan"%in%sc$path){ #if in pan
    print("sp at pan")
    selection <- pandata[pandata$X.1%in%spi,]
    t <- colSums(matrix(as.logical(as.numeric(as.matrix(selection[,-1]))), nrow=nrow(selection)))
    #set all NA's to 1
    t[is.na(t)] <- 1     
    if (length(t)!=0&sum(t)>0){
      t <- as.logical(t)
      # occzones <- subset(zones, t)
      occzones <- zones[zones@data$PanArticZo%in%names(selection)[-1][t],]
      #see if there are overlaps with hultenGBIF data
      temp_rasta <- rasta
      temp_rasta[][is.na(temp_rasta[])] <- 0 #set NA's to zero
      upa <- extract(temp_rasta, occzones)
      upa <- unlist(lapply(upa, max, na.rm=T))
      if (any(!upa)){ # if we endup nowhere,,,
        occzones2 <- subset(occzones, !upa)
        # create a raster from occzones...
        occzonesrasta <- rasterize(occzones2, disaggregate(rasta, fact=10*1), field=1)
        occzonesrasta <- aggregate(occzonesrasta, fact=10*1)
        # merge final
        rasta <- merge(occzonesrasta, rasta)
        ## plot(crop(occzonesrasta, extent(-180, 180, 10, 90)), legend=FALSE, axes=FALSE, box=FALSE, add=T, col=rgb(0,0,1,0.5))
      }
    }
  }
  #plot
  ## plot(crop(rasta, extent(-180, 180, 10, 90)), main=paste("FINAL: ", spi), legend=FALSE, axes=FALSE, box=FALSE)
  ## dev.off()
  savetodisk(raster=rasta, sp_name=spi, output_dir=output_dir_res)
}
