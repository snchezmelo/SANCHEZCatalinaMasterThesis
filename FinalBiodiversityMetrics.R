####This script allows to calculate biodiversity metrics of SR, PD, PE from SDMs.
## Original from Andrea PAZ, adapted by Catalina SANCHEZ.

##################
#Species richness#
##################


library(raster)
library(dismo)

#Load binary reclassified raster files for individual species
models <- list.files(path = "/Volumes/TOSHIBA_EXT/RecclasifiedPass/", pattern=".tif")

raster_list <- lapply(models, raster)

rasters_stack <- stack(models)

#Calculate species richness by summing across the stack
species_richness <- calc(rasters_stack, sum)

#Plot the species richness map
plot(species_richness, main = "Species Richness")

#Save the species richness raster
writeRaster(species_richness, filename = "species_richness_AllPass.tif", format = "GTiff", overwrite = TRUE)

#########################
#Phylogenetic Diversity#
#########################

library(rgdal)
library(raster)
library(picante)
library(geiger)


models<-list.files(path = "/Volumes/TOSHIBA_EXT/ReclassifiedPassForPD/Resampled10KmReclassifiedPassForPD/", pattern=".tif")

#Process the species names without the .tif extension
species_names_full <- sub(".tif$", "", basename(models))


####Load Phylogeny, generate list of species
##To read a phylogeny in newick format use read.tree instead of read.nexus
#In working directory always a file with the same name: phylogeny.nex if not then change name before

##read phylogeny
user_phylogeny<-ape::read.tree("KozakAdjusted3.tree")

####MMake sure the names match between maps and phylogeny
setdiff(user_phylogeny$tip.label, species_names_full)
##if it doesnt use this code
###Trim phylogeny to match distribution data (remove non hylids)
pruned.tree<-ape::drop.tip(user_phylogeny, setdiff(user_phylogeny$tip.label, species_names_full))
##check that all disrtribution data is in tree
test<-as.data.frame(species_names_full)
rownames(test)<-species_names_full
check_names<-geiger::name.check(user_phylogeny, test, data.names=NULL) #this must be OK
### Trim the distribution data to match phylogeny
#species_names_full<-species_names_full[species_names_full %in% user_phylogeny$tip,label]
#distribution_files<-sub("*$","_OR_AUC_T10.asc",species_names)

## Now we can actually start the analysis

library(terra)

tiff_files <- list.files(path = "/Volumes/TOSHIBA_EXT/ReclassifiedPassForPD/Resampled10KmReclassifiedPassForPD/", pattern=".tif")

species_names_full <- user_phylogeny$tip.label
Stack_full <- rast(tiff_files)

tabla<-as.data.frame(species_names_full)
colnames(tabla)<-"Grilla"

#Create empty raster of the desired area and resolution to assign pixel numbers
r<-Stack_full[[1]]
res(r)<-res(Stack_full) #resolution
grilla=r
names(grilla)="Grilla"
grilla[1:ncell(grilla)]<-1:ncell(grilla)
Stack <- c(grilla,Stack_full)

#Turn maps into dataframe for computation of PD
community_data<-terra::as.data.frame(Stack)
##remove rows that are NA 
#The 171 number should change depending on the amount of models, here it was 170
community_data1<-community_data[-which(rowSums(!is.na(community_data[,2:171]))==0),]
species_names<-colnames(community_data1)[2:length(community_data1)] #Store species names

#Store 
write.table(community_data1, file = "communities_PassifloraForPD.txt", append = FALSE,row.names=F,quote=F,sep="\t") #tratar de agregar nombre de mascara

#Phylogenetic diversity computation 

pd.result <-picante::pd(community_data1[,2:ncol(community_data1)],user_phylogeny,include.root = F) 
#Add the pixel PD value to data frame
community_data1$pd<-pd.result[[1]]
#And we also add SR:
community_data2 <-  community_data1
community_data2$sr<-pd.result[[2]]
#Write the new matrix to a file to avoid rerunning the script for potential further analyses
write.table(community_data2, file = "communities_PassifloraForPD_with_pd_and_sr.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#Generate a raster containing PD information per pixel

#1-First generate an empty raster using a base model for resolution and area

values(r)<-0
pd_ras<-r
values(pd_ras)<-NA #Eliminate every value they will be replaced by PD values further down


#2- Assign PD values to raster
pd_ras[community_data1$Grilla]<-community_data1$pd

#3- Save raster to file 

terra::writeRaster(pd_ras,"Passiflora_PD.asc")

#4-Optional plotting map in R 
plot(pd_ras)

#Phylogenetic endemism

pe_values <- picante::pd(community_data1, user_phylogeny)$PD / rowSums(community_data1)  # This is a simplified calculation
pd_values <- pd.result$PD
# Calculate PE as PD divided by the number of species per grid cell (row)
pe_values <- pd_values / rowSums(community_data_filtered)

#Create an empty raster based on your template
template_raster <- raster_stack[[1]]
pe_raster <- template_raster
values(pe_raster) <- NA  # Initialize with NA

#Assuming your community data has a 'grilla' column with pixel indices
pe_raster[community_data1$Grilla] <- pe_values

#Save the PE raster to a file
writeRaster(pe_raster, filename = "PE_resultPass.tif", format = "GTiff", overwrite = TRUE)

#Optional: Plot the PE raster
plot(pe_raster)
