# Species_Range_Mapping

## **Title** 
Species Range Mapping
R-scripts to generate species range from occurance data.

## **Author**
Oskar Hagen (oskar@hagen.bio)
for range map function: Oskar Hagen, Fabian Fopp and Camille Albouy

## **Description** 
In a simplified way, the four main process from data to species ranges are: 
1.	Extract occurrence data from GBIF species that are pre-selected as cold-adapted or that complement Hulten’s ranges maps and the Panarctic Flora checklist. For that, data manipulation and a pre-resolution of species taxonomic names is necessary (see ‘1.extract_pre-select_cold_plants_gbif.R’).
2.	After resolving taxonomic names and synonyms within datasets, merge distribution data (see ‘2.merge_dataset.R’).
3.	From all merged points, create range maps using ecoregions and climatic layers information (see ‘3.create_range.R’).
4.	After having clean datasets, merge all datasets (see ‘4.merge_datasets.R’).
Additional steps as the removal of plants not belonging to the Angiosperms clade, masking of the cold areas and the manual check of species names are not provided here, given their simplicity, specificity to a cluster computer and large quantity of intermediate individual files.
Inside the folder ‘support_functions’, there are four functions called by the scripts mentioned above. All functions have a short description and parameters are shortly described.
The file (i) ‘range_from_points.R’ contains the function used to create range maps from occurrence data. The file (ii) ‘resolve_taxonomic_names.R’ contains the function used to automatically resolve taxonomic names and look for synonyms while deleting clear mistakes. The taxonomic name resolution procedure involved a manual check of all species. The file (iii) ‘merge_occurance_data.R’ merges occurrence data from a same dataset considering duplicates and corrected species names. This function requires a species names resolution matrix and a specific structure for datasets (i.e. occurrences inside a folder named ‘points’ with a txt file for each species name containing x and y coordinates). The final occurrence data is stored in a folder named ‘merged_point’. The file (iv) ‘compress_housekeeping.R’ contains a function conveniently called after generating the merged points. It compresses all occurrence data inside ‘point’ folder and deletes the folder while keeping a compressed zip file for reference. 

## **Reference** 
Oskar Hagen, Lisa Vaterlaus, Camille Albouy, Andrew Brown, Flurin Leugger, Renske E. Onstein, Charles Novaes de Santana, Christopher R. Scotese, Loïc Pellissier. (2019) Mountain building, climate cooling and the richness of cold-adapted plants in the Northern Hemisphere. Journal of Biogeography
