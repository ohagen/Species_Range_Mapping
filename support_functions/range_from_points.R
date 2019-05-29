##################################################
## Description: Create range maps from points based on 
## eco/bio-regions, climate and outliers
## 
## Date: 2018-05-21 13:26
## Authors: Oskar Hagen (oskar@hagen.bio)
##         Fabian Fopp (fabian.fopp@usys.ethz.ch)  
##         Camille Albouy (albouycamille@gmail.com)
##################################################

range.fun <- function (species_name, occ_coord, proj, Climatic_layer, Bioreg, final_resolution, degrees_outlier=3, clustered_points_outlier=1,
                       buffer_width_point=0.5, buffer_increment_point_line=0.5, buffer_width_polygon=0.1, cut_off=0.05, method="PC", cover_threshold=0.3,
                       write_plot=F, dest_plot='./Output/Maps/', write_raster=T, dest_raster='./Output/Rasters/', return_raster=T, overwrite=F, dir_temp='.', dir_log='./log/'){
  
  #Function description
  # This function estimates species ranges based on occurrence data, eco/bio-regions and a climatic layer. It first deletes outliers from the 
  # observation dataset. Then creates a polygon (convex hull) with a user specified buffer around all the observations of one bioregion. 
  # If there is only one observation in a bioregion, a buffer around this point is created. If all points in a bioregion are on 
  # a line, the function will also create a buffer around these points but the buffer size increases with the number of points in the line
  
  
  #Parameters description
  # species_name: character string with the species name. E.g. "Anemone nemorosa".
  # occ_coord: a dataframe with two columns containing the coordinates of all observations of a given species. 
  # proj: Spatial projection in which the coordinates of the occurrence data (input) are stored. The output raster of the species range will have the same projection.
  # Climatic_layer: climate raster (e.g. temperature) used to improve the distribution range (by rejecting cells with unsuitable climate).
  # Bioreg: shapefile containg different eco/bio-regions that will be used for defining convex hulls.
  # final_resolution: determines the final resolution of the species range raster.
  # degrees_outlier: distance threshold (degrees) for outlier classification. If the nearest minimal distance to the next point is larger than this threshold, it will be considered as an outlier. 
  # clustered_points_outlier: maximum number of grouped points (closer than the degrees_outlier to each other) that should still be considered as outliers if isolated from all other points (given degrees_outlier).
  # buffer_width_point: buffer (in degrees) which will be applied around single observations.
  # buffer_increment_point_line: how much should the buffer be increased for each point on a line given in percentage.
  # buffer_width_polygon: buffer (in degrees) which will be applied around distribution polygons (for each bioregion)
  # cut_off: quantile of temperatures which are not considered to be suitable for the species. e.g: cut_off=0.05: lowest and highest 5% of temperature distribution
  # method: method used to define presences and absences, can be: 'PC' (percentage cells method): uses temperature filter only at the edge of range raster or 'AC' (all cells method) method for all cells. See methods details for more information
  # cover_threshold: only relevant if method="PC". Specifies the proportion threshold of covered area by the range polygon above which the corresponding raster will be classified as "present" without considering the temperature.
  # write_plot: should a range plot saved? (TRUE or FALSE)
  # dest_plot: path to where the plot should be saved (only relevant if write_plot=TRUE)
  # write_raster: should the range raster be saved? (TRUE or FALSE)
  # dest_raster: path to where the raster should be saved (only relevant if write_raster=TRUE)
  # return_raster: should the raster be returned by the function? (TRUE or FALSE)
  # overwrite: if the plot/raster for this species already exists, should it be overwritten? (TRUE or FALSE)
  # dir_temp: location for temporary storage of convex hull text file analysis?
  # dir_log: location of log files (e.g. to report insufficient occurrences etc.)?
  
  
  ##### METHODS DETAILS #####
  #PC (percentage cells method)
  # In the conversion from the spatial polygons to the final range raster, every raster cell that is covered by a certain threshold (default=0.3)
  # will be marked as present. If the cell is covered by less than this threshold (but higher than 0), the cell will be marked as absence
  # if the climatic value of this cell is below or above the climatic quantile declared. 
  #Otherwise the cell will be marked as present.
  #AC (all cells method)
  # In the conversion from the spatial polygons to the final range raster, the cell will be marked as absence if the climatic value of 
  # this cell is below the 0.05 or above the 0.95 climatic quantile of all observations (quantile can be changed). 
  # Otherwise the cell will be marked as present.
  
  
  # function relevant packages
  lib_vect <- c("raster","rgdal","sp","maptools","rgbif","shape","geometry","rgeos", "FNN", 'maps')
  sapply(lib_vect,require,character.only=TRUE)
  
  # handle climatic layer
  if(class(Climatic_layer)=='RasterLayer'){
    if(final_resolution<res(Climatic_layer)[1]){
      warning('Resolution of climatic layer (', round(res(Climatic_layer)[1],3), ') is coarser than desidered final resolution (', final_resolution, ').
              ', round(res(Climatic_layer)[1],3), ' will be used as final resolution')
      final_resolution <- res(Climatic_layer)[1]
    }
  }
  
  
  # grid making
  grd <-  SpatialGrid(GridTopology(cellcentre.offset=c( -179.5, -89.5), cellsize=c(final_resolution, final_resolution),cells.dim=c(360,180)/final_resolution), proj4string=proj)
  
  
  if(class(Climatic_layer)=='RasterLayer'){
    Climatic_layer <- projectRaster(Climatic_layer, raster(grd)) 
  }
  
  #remove duplicates. #here depending on desired precision, rounding is recomended to avoid future problems with the convex hull analysis
  occ_coord <- unique(occ_coord) 
  
  
  #### START INTERNAL FUNCTIONS DECLARATION ####
  conv_function <- function (x=coord_2,polygons=TRUE,proj=proj){
    
    if (nrow(x)<3){ #check number of observations points in each bioregion, if <3 create point buffer
      sp_coord <- SpatialPoints(x,proj4string=proj)
      return(gBuffer(sp_coord,width=buffer_width_point+(nrow(x)-1)*buffer_increment_point_line)@polygons[[1]]) 
      
    } else {
      #test if points are on line, if yes create point buffer around points
      is_line <- 0
      for(i in 2:(nrow(x)-1)){
        # y <- round(x,3) #to avoid problems with colinearity in polygon (like that it will be classified as a line)
        dxc <- x[i,1]-x[i-1,1]
        dyc <- x[i,2]-x[i-1,2]
        dx1 <- x[i+1,1]-x[i-1,1]
        dy1 <- x[i+1,2]-x[i-1,2]  
        is_line[i-1] <- dxc*dy1-dyc*dx1
      }
      
      if(all(abs(is_line)==0)){ #or better round here?
        cat('Bioreg=', g, nrow(x), 'points laying on one line. Using buffer width of ', buffer_width_point+(nrow(x)-1)*buffer_increment_point_line, '\n')
        sp_coord <- SpatialPoints(x,proj4string=proj)
        return(gBuffer(sp_coord,width=buffer_width_point+(nrow(x)-1)*buffer_increment_point_line)@polygons[[1]])
      }
      else { #if they are not on a line, create a convex hull
        arg_file <- paste0('QJ Fx TO ', file.path(dir_temp, 'vert.txt'))
        vert0<-convhulln(x, arg_file)
        vert1<-scan(file.path(dir_temp, 'vert.txt'),quiet=T);file.remove(file.path(dir_temp, 'vert.txt'))
        vert2<-(vert1+1)[-1]
        FE_vert<-row.names(x)[vert2]
        coord_conv <- x[FE_vert,]
        
        if(polygons==TRUE){
          coord_conv <-rbind(coord_conv,coord_conv[1,])
          P1 <- Polygons(srl=list(Polygon(coord_conv,hole=FALSE)),ID="PolygA")
          P1 <- SpatialPolygons(Srl=list(P1),proj4string=proj)
          return(gBuffer(P1,width=buffer_width_polygon)@polygons[[1]]) #CA
        } else {
          return(coord_conv)
        } # end of ifelse
      } #end of ifelse  
    }# end of ifeslse
  } # end of conv_function
  
  plot.occ <- function(){
    if(!dir.exists(paste(dir_log, 'plots', sep='/'))){
      dir.create(paste(dir_log, 'plots', sep='/'), recursive=T)
    }
    jpeg(paste0(dir_log, 'plots/', species_name, '.jpg'), width=3000, height=2000)
    plot(Bioreg)
    points(occ_coord, col='red', pch=4, cex=1.5, lwd=1.2)
    points(occ_coord, col='red', pch=0, cex=1.5, lwd=1.2)
    dev.off()
  }
  
  #### END INTERNAL FUNCTIONS DECLARATION ####
  
  if (nrow(occ_coord)<=clustered_points_outlier+1){
    plot.occ()
    warning('Too few occurrences.')
    writeLines(c('Too few occurrences', species_name), con=paste0(dir_log, species_name, '.txt'), sep=' ')
  } else{ 
    
    cat("########### Start of computation for species: ",species_name," #######################", "\n") 
    
    #create distance matrix...
    mat_dist <- as.matrix(knn.dist(occ_coord, k=clustered_points_outlier))
    
    
    #mark outliers
    cond <- apply(mat_dist, 1, function(x) x[clustered_points_outlier])>degrees_outlier
    rm(mat_dist) 
    
    print(paste0(sum(cond), " outlier's from " ,nrow(occ_coord), " | proportion from total points: ", round((sum(cond)/nrow(occ_coord))*100,0), "%"))
    
    occ_coord_mod <- occ_coord[!cond,]
    
    if(nrow(occ_coord_mod)==0){
      warning('Too few occurrences within outlier threshold.')
      writeLines(c('Too few occurrences within outlier threshold for', species_name), con=paste0(dir_log, species_name, '.txt'), sep=' ')
      plot.occ()
      cat("########### End of computation for species: ",species_name," #######################", "\n") 
      return(NULL)
      
    } else{
      
      #handy debug plot line
      #map('world'); points(occ_coord, col='red', pch=15); points(occ_coord_mod, col='blue', pch=16)
      
      #correct rownames 
      rownames(occ_coord_mod) <- 1:dim(occ_coord_mod)[1]
      
      occ_points <- SpatialPoints(coords=occ_coord_mod,proj4string=proj)
      
      cat("### Projection adjustement for bioregion shapefile...", "\n") 
      unique <- unique(Bioreg@data$ECO_NAME)
      gc()
      
      cat("### Interscetion between occurences and bioregions ...", "\n") 
      
      SP_dist <- list()
      a <- 0
      for(g in 1:length(unique)) {
        tmp <- as(gSimplify(Bioreg[Bioreg$ECO_NAME == unique[g],],tol=0.001,topologyPreserve=TRUE),"SpatialPolygons")
        a <- data.frame(occ_coord_mod[names(na.omit(over(occ_points,tmp,fn=NULL))),, drop=F])
        if (nrow(a)==0) {
          SP_dist[[g]] <- NA
        } else {
          #cat("g=",g,'\n')
          SP_dist[[g]] <- suppressWarnings(gIntersection(gBuffer(SpatialPolygons(Srl=list(conv_function(a,proj=proj))), width=0),tmp)) #zero buffer to avoid error
          if(class(SP_dist[[g]])=='SpatialCollections'){
            SP_dist[[g]] <- SP_dist[[g]]@polyobj #only keep SatialPolygon
          }
          
        } # end of if
        
      } # end of SP_dist
      
      L <- SP_dist[!is.na(SP_dist)]
      
      if(length(L)==0){
        plot.occ()
        warning('No occurrences within Bioregions. Empty raster produced.')
        writeLines(c('No occurrences within Bioregions. Empty raster produced for', species_name), con=paste0(dir_log, species_name, '.txt'), sep=' ')
        return(NULL)
      } else{
        
        names_poly <- paste("Poly_",seq(1,length(L),1),sep="")
        for(i in 1:length(L)){L[[i]]@polygons[[1]]@ID <- names_poly[i]}
        
        shp_species <- SpatialPolygons(Srl=lapply(L, function(x){x@polygons[[1]]}),proj4string = proj)
        
        if (method=='PC'){ 
          cat("### Using Percentage Cells Method ###", "\n")
          range_raster <- rasterize(shp_species, raster(grd), getCover=T)
          range_dataframe <- as.data.frame(range_raster, xy=T)
          
          if(class(Climatic_layer)=='RasterLayer'){
          
          
          temp_all <- extract(Climatic_layer, occ_coord_mod[, 1:2]) #store all temperatures at observed points
          moderate_temperature <- temp_all[temp_all<quantile(temp_all,probs=1-cut_off,na.rm=TRUE) & 
                                             temp_all>quantile(temp_all,probs=cut_off,na.rm=TRUE)] #store moderate temperature
          moderate_temperature <- moderate_temperature[!is.na(moderate_temperature)]
          range_dataframe$temp <- extract(Climatic_layer, range_dataframe[, c('x', 'y')]) #extract temperature of all raster cells
          
          range_dataframe$occ <- ifelse(range_dataframe$layer>=(cover_threshold * 100), 1, 
                                        ifelse(range_dataframe$temp>=min(moderate_temperature) & range_dataframe$temp<=max(moderate_temperature) & range_dataframe$layer>0, 1,
                                               ifelse(is.na(range_dataframe$temp), NA, 0)))
          
          
          } else { 
            range_dataframe$occ <- ifelse(range_dataframe$layer>=(cover_threshold * 100), 1, 0)
            }
          rast_end <- rasterFromXYZ(range_dataframe[, c('x', 'y', 'occ')], res=final_resolution)
        } 
        else if (method=='AC') { 
          cat("### Using All Cells Method ###", "\n")
          pts_d <- do.call(rbind,lapply(shp_species@polygons,function(z){z@Polygons[[1]]@coords}))
          
          if(class(Climatic_layer)=='RasterLayer'){
            a <- extract(Climatic_layer,occ_points)
            
            test <- Climatic_layer>quantile(a,probs=cut_off,na.rm=TRUE)
            test_2 <-  Climatic_layer<quantile(a,probs=1-cut_off,na.rm=TRUE ) 
            yop <- test_2+test
          }
          
          cat("### Raster making ###", "\n")
          grd <- raster(grd); grd[] <-0
          rast_cel <- unique(cellFromXY(grd,pts_d))
          data <- do.call(rbind,cellFromPolygon(grd,shp_species,weights=TRUE))[,1]
          a<-grd; a[] <- 0
          a[][unique(c(data,rast_cel))]<-1
          
          if(class(Climatic_layer)=='RasterLayer'){
          yop_2 <- crop(yop,a) 
          yop_2 <- yop_2>1
          
          rast_end <- ((yop_2 +a)==2)
          } else{
            rast_end <- a
          }
          
          
        } else{
          stop("Please choose a valid method!")
        }
      
        #plot species range
        if(write_plot==T){
        if(!dir.exists(dest_plot)){
          dir.create(dest_plot, recursive=T)
        }
        if(file.exists(paste0(dest_plot, species_name, '.jpg')) & overwrite==F){
          stop('Plot already exists. Use overwrite=T to overwrite it.')
        }
        jpeg(paste0(dest_plot, species_name, '.jpg'), height=2000, width=3000)
        plot(rast_end, col=c(rgb(0,0,0,0), rgb(0,1,0,1)), legend=F) 
        map('world', add=T)
        dev.off()
      }
        
      if(write_raster==T){
        if(!dir.exists(dest_raster)){
          dir.create(dest_raster, recursive=T)
        }

        if(file.exists(file.path(dest_raster, paste0(species_name, ".asc"))) & overwrite==F){
          stop('Raster already exists. Use overwrite=T to overwrite it.')
        }
        writeRaster(rast_end,filename=file.path(dest_raster, paste0(species_name, ".asc")), overwrite=TRUE)
      }
        
        if(return_raster==T){
          cat("########### End of computation for species: ",species_name," #######################", "\n") 
          return(rast_end)
        }
        rm(rast_end)
      } #end else L!=list()
    }# ende else 'enough occurrences'
  }
}