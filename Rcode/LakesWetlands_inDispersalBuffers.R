############### Lake and wetland patches within animal dispersal buffers #######################
# Date: 12-15-18
# updated: 4-22-19
# Author: Ian McCullough, immccull@gmail.com
################################################################################################

#### R libraries ####
library(raster)
library(rgeos)
library(spatialEco)

#### Input data ####
setwd("C:/Users/FWL/Documents/LivinOnTheEdge")

# Michigan shapefile
mich_shp <- shapefile("Data/GIS/Michigan_NoIsleRoyale.shp")

# Michigan lakes from NHD (and convert to same crs as other data) #uncut borders allows lakes to spill across state line instead of being cut off (selected spatially in ArcGIS)
# USGS. (2018). National Hydrography Dataset Plus medium resolution version 2. Accessed 
# from: https://www.epa.gov/waterdata/get-data in January 2019
mich_lakes_NHD <- shapefile("C:/Ian_GIS/NHD/NHD_waterbody_Mich/NHD_waterbody_Mich_unCut_borders.shp") 
mich_lakes_NHD <- spTransform(mich_lakes_NHD, CRSobj=crs(mich_shp))

# LAGOS GIS data downloaded and stored locally from (including mich_shp above): 
# Soranno P., K. Cheruvelil. (2017). LAGOS-NE-GIS v1.0: A module for LAGOS-NE, 
# a multi-scaled geospatial and temporal database of lake ecological context and water 
# quality for thousands of U.S. Lakes: 2013-1925. Environmental Data Initiative. 
# Package ID: edi.98.1
# http://dx.doi.org/10.6073/pasta/fb4f5687339bec467ce0ed1ea0b5f0ca. Dataset accessed 9/26/2017.
# LAGOS NE lakes
#lakes_4ha_pts <- shapefile("C:/Ian_GIS/LAGOS_US_4ha_lakes/LAGOS_US_All_Lakes_4ha_pts/LAGOS_US_All_Lakes_4ha_pts.shp")
lakes_4ha_poly <- shapefile("C:/Ian_GIS/LAGOS-NE-GISv1.0/LAGOS_NE_All_Lakes_4ha/LAGOS_NE_All_Lakes_4ha.shp")

# NWI polygons (selected in ArcGIS from MI,WI,OH and IN; polygons within 5000m of Mich border to allow for connectivity outside of Michigan)
# WETLAND_TYPE = 'Freshwater Emergent Wetland' OR WETLAND_TYPE = 'Freshwater Forested/Shrub Wetland' OR WETLAND_TYPE = 'Other'
# USFWS. (2018). National Wetlands Inventory website. U.S. Department of the Interior, Fish and 
# Wildlife Service, Washington, D.C. http://www.fws.gov/wetlands/ (accessed February 2019)
mich_NWI <- shapefile("C:/Ian_GIS/FreshwaterConservation/NWI/MI_NWI_5000mBuff/MI_NWI_5000mBuff.shp")

#### D-fine function ####
# Calculates statistics for wetlands and small lakes within dispersal buffer around focal lake
source("Rcode/functions/lake_dispersal_buffer_patch_metrics.R")

####################### Main program #############################
# identify Mich lagoslakeids (focal lakes)
mich_lakes_4ha <- subset(lakes_4ha_poly, STATE=='MI')
mich_lagoslakeids <- unique(mich_lakes_4ha@data$lagoslakei)

# get rid of Great Lakes from NHD
great_lakes_comid <- c(904140245, 904140244, 904140243, 904140248)
mich_lakes_NHD_sub <- subset(mich_lakes_NHD, !(COMID %in% great_lakes_comid))

# create empty data frame for loop, with each iteration filling in a row
LakeBufferPatchStatz <- data.frame(matrix(NA, nrow = length(mich_lagoslakeids), ncol = 10))
colnames(LakeBufferPatchStatz) <- c('lagoslakeid','nLakePatches','LakeEdge_km','LakeEdgeArea_ha','LakeEdgeArea_pct',
                                    'nWetlandPatches','WetlandEdge_km','WetlandArea_ha','WetlandArea_pct','BufferWidth_m')
dispersal_buff <- 2020 #meters from Patrick et al. (2012) for common snapping turtle
#Patrick, D. A., Gibbs, J. P., Popescu, V. D., & Nelson, D. A. (2012). Multi-scale habitat-resistance models for predicting road mortality "hotspots" for turtles and amphibians. 
#Herpetological Conservation and Biology, 7(3), 407-426.

# loop took almost 2 days fro 6511 Michigan lakes
for (i in 1:length(mich_lagoslakeids)){
  tump <- lake_dispersal_buffer_patch_metrics(lagoslakeid=mich_lagoslakeids[i], LAGOS_shp=lakes_4ha_poly, NHD_shp=mich_lakes_NHD_sub, wetland_shp=mich_NWI, dispersal_buff=dispersal_buff)
  LakeBufferPatchStatz[i,] <- tump
  tump <- NULL
}

# correct negative lake metric values that arose from subtracting focal lake perimeter (convert negatives to 0)
LakeBufferPatchStatz$LakeEdge_km <- ifelse(LakeBufferPatchStatz$LakeEdge_km < 0, 0, LakeBufferPatchStatz$LakeEdge_km)
LakeBufferPatchStatz$LakeEdgeArea_ha <- ifelse(LakeBufferPatchStatz$LakeEdgeArea_ha < 0, 0, LakeBufferPatchStatz$LakeEdgeArea_ha)
LakeBufferPatchStatz$LakeEdgeArea_pct <- ifelse(LakeBufferPatchStatz$LakeEdgeArea_pct < 0, 0, LakeBufferPatchStatz$LakeEdgeArea_pct)

#write.csv(LakeBufferPatchStatz, file="Data/LakeWetlandPatchStats_2020mBuff.csv")

par(mfrow=c(2,3))
hist(LakeBufferPatchStatz$nLakePatches)
hist(LakeBufferPatchStatz$LakeEdge_km)
hist(LakeBufferPatchStatz$LakeEdgeArea_ha)
hist(LakeBufferPatchStatz$LakeEdgeArea_pct)
hist(LakeBufferPatchStatz$nWetlandPatches)
hist(LakeBufferPatchStatz$WetlandEdge_km)
hist(LakeBufferPatchStatz$WetlandArea_ha)
hist(LakeBufferPatchStatz$WetlandArea_pct)
