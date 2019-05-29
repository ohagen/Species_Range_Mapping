##################################################
## Description: Select and extract cold species
## from gbif downloaded data:
## GBIF.org (26th August 2018) GBIF Occurrence Download https://doi.org/10.15468/dl.1xd4qs
## As well as overlapping species with Hulten's and PAF are also extracted if present
## since we consider all species in PAF and Hulten's cold adapted at this initial stage
##
## Date: 2018-09-19 21:22:02
## Author: Oskar Hagen (oskar@hagen.bio)
##################################################

#load libraries
library(data.table)
library(raster)
library(Taxonstand)

#define variables
time <- "20181024"
new.extent <- extent(-180, 180, 10, 90)

#set working directory
setwd("SET YOUR WORKING DIRECTORY HERE")

#load temperature map HIGH resolution i.e. 2.5 minutes (avaiable at Worldclim)
climate<-raster("YOUR CLIMATE LAYER HERE")
#temperature for north hemisphere wusing extent format (xmin,xmax,ymin,ymax)
climate <- crop(x = climate, y = new.extent)

#in order to handle this huge amount of data, a first selection of tables was done using a UNIX machine
#below shell code to select big data cols:

#shell#   awk -F'\t' '{print $10,"\t",$17,"\t",$18}' 0000680-180824113759888.csv > subset_tab.csv
   
#R#       t <- data.table:::fread(input="G:/data/raw_gbif/subset_tab.csv")
#R#       t <- t[species!=""]
#R#       t <- t[species!=""]
#R#       #write
#R#       fwrite(t, file="G:/data/raw_gbif/subset_tab_clean.csv")

#load
t <- data.table:::fread(input="C:/VITAL LOCAL/Meus Documentos/ETH PhD/MainPapers/ColdPlants/data/GBIF_tracheophytes/subset_tab_clean.csv")
#this will do the job and get the 50 quantile temperature.
t[, temperature := quantile(round(extract(climate, .SD ),2), 0.75, na.rm=T), by=.(species), .SDcols = c("decimallongitude", "decimallatitude")]
#get cold species
coldsp <- unique(t[temperature<=5]$species)

#we now make sure we get the points of species that are in hulten or paf
#get ALL gbif species ---------------
spgbif <- unique(t$species)
#tpl it for reference
tgbif <- TPL(spgbif, diffchar = 4, drop.lower.level=T)
mgbif <- cbind(submittedname=spgbif, finalname=paste(tgbif[,"New.Genus"], tgbif[,"New.Species"]))
#save
saveRDS(mgbif, file="./MainPapers/ColdPlants/scripts/analysis/resolution_1_gbif.RDS")
mgbif <- readRDS(file="./MainPapers/ColdPlants/scripts/analysis/resolution_1_gbif.RDS")



#prepare list
sdf <- list()

#get Hulten data ----------------------
h_dir <- "./MainPapers/ColdPlants/RetriveDistData/2_hulten/point"
h_files <- list.files(h_dir)
length(h_files)
h_cold_sps <- NULL
for (i in 1:length(h_files)){
  p <- read.table(file.path(h_dir, h_files[i]), header = T)
  q <- quantile(round(extract(climate, p ),2), 0.50, na.rm=T)
  print(q)
  if (q<=5.1&!is.na(q)){
    print(paste("this species", h_files[i], "was below 5.1, by: ", q))
    h_cold_sps <- c(h_cold_sps, tools::file_path_sans_ext(h_files[i]))
  }
}
h_all_sps <- gsub(".txt$", "", h_files)
#get list of species to discart...
h_noncold_sps <- h_all_sps[!h_all_sps%in%h_cold_sps]

sdf[[1]] <- cbind("2_hulten", h_cold_sps)

#get paf data ----------------------
pandata <- readRDS("./MainPapers/ColdPlants/scripts/RetriveDistData/4_panarticflora/scripts/presences(PANARTIC).rds")
paf_species <- pandata$X.1
removethis <- c(" x$", " \"$", " sp.$", " ssp.$", "\\d", "\\b\\. ", "^NA ", "^[a-zA-Z]{2} ","^[a-zA-Z]{1} ", "^NA$", " NA$", "'", "\"", "/")
# mistakes <- grep(paste(removethis,collapse="|"), paf_species, value=T)
# grep(" x$", paf_species, value=T, invert=T)
paf_species <- grep(paste(removethis,collapse="|"), paf_species, value=T, invert=T)

write.csv2(paf_species, file="./MainPapers/ColdPlants/scripts/analysis/PAF_species_pre_filtered.csv", row.names = F)

sdf[[2]] <- cbind("4_panarticflora", paf_species)

rm(pandata)

#merge hulten and paf --------------------
sdf <- do.call(rbind, sdf)
colnames(sdf) <- c("path", "submittedname")
sphp <- unique(sdf[,"submittedname"])
#tpl it all for reference
thp <- TPL(sphp , diffchar = 4, drop.lower.level=T)
mhp <- cbind(submittedname=sphp, finalname=paste(thp[,"New.Genus"], thp[,"New.Species"]))
elements <-   unlist(lapply(strsplit(mhp[, "finalname"], "\\s+"), length))
mhp <- mhp[!elements!=2, ,drop=F]
#remove NA's
mhp <- na.omit(mhp)
#save
saveRDS(mhp, file="./MainPapers/ColdPlants/scripts/analysis/resolution_3_4_hulten_paf.RDS")

mhp <- readRDS(file="./MainPapers/ColdPlants/scripts/analysis/resolution_3_4_hulten_paf.RDS")



### NEW REMOVAL OF NON COLD GIVEN hulten and PAF ####

#to remove hulten
h_toremove_sps <- mhp[mhp[,1]%in%h_noncold_sps,2]
length(h_toremove_sps)
h_toremove_sps <- unique(h_toremove_sps)

#to remove PAF
final_gbif_warm <- mgbif[!mgbif[,1]%in%coldsp,2]
warm_paf <- mhp[mhp[,2]%in%final_gbif_warm,2]
paf_toremove_sps <- warm_paf

to_remove_h_paf <- c(h_toremove_sps, paf_toremove_sps)
to_remove_h_paf <- unique(to_remove_h_paf)

#save
saveRDS(h_toremove_sps, file=paste0("./MainPapers/ColdPlants/scripts/analysis/", time,"_h_toremove_sps.RDS"))
saveRDS(paf_toremove_sps, file=paste0("./MainPapers/ColdPlants/scripts/analysis/", time,"_paf_toremove_sps.RDS"))
saveRDS(to_remove_h_paf, file=paste0("./MainPapers/ColdPlants/scripts/analysis/", time,"_to_remove_h_paf.RDS"))


# overlap missed -----------------
#get overlap
om <- mgbif[mgbif[,"finalname"]%in%mhp[,"finalname"],]

#remove already considered cold by gbif only
additional <- om[!om[,"submittedname"]%in%coldsp, ]

extractgbif <- unique(c(coldsp, additional[,"submittedname"]))

#subset cold plants...
pb <- txtProgressBar(min = 0, max = length(extractgbif), style = 3)
i <- 1
for (sp_i in extractgbif){
  #print(paste("working with: ", sp_i))
  sp_occ <- t[grepl(sp_i, species)][,c(3,2)] #slower but precise in case of problems with on species names...
  # sp_occ <- t[species==sp_i][,c(3,2)] #faster but requires the treatment of subspecies later on
  names(sp_occ) <- c("X", "Y")
  #print("saving...")
  # write.table(sp_occ, file = paste0("C:/VITAL LOCAL/Meus Documentos/ETH PhD/MainPapers/ColdPlants/data/GBIF_tracheophytes/block/", sp_i, ".txt"), row.names = FALSE)
  write.table(sp_occ, file = paste0("./MainPapers/ColdPlants/data/GBIF_tracheophytes/block/bypass_filter_2/", sp_i, ".txt"), row.names = FALSE)
  i <- i+1
  setTxtProgressBar(pb, i)
}
close(pb)