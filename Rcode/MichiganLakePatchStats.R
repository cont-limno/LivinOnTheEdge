####################### Patch stats on Michigan lakes ##########################################
# Date: 11-8-18
# updated: 4-22-19
# Author: Ian McCullough, immccull@gmail.com
################################################################################################

#### R libraries ####
library(raster)
library(SDMTools)

###### Input data ######
setwd("C:/Users/FWL/Documents/LivinOnTheEdge")

# LAGOS GIS data downloaded and stored locally from (including mich_shp above): 
# Soranno P., K. Cheruvelil. (2017). LAGOS-NE-GIS v1.0: A module for LAGOS-NE, 
# a multi-scaled geospatial and temporal database of lake ecological context and water 
# quality for thousands of U.S. Lakes: 2013-1925. Environmental Data Initiative. 
# Package ID: edi.98.1
# http://dx.doi.org/10.6073/pasta/fb4f5687339bec467ce0ed1ea0b5f0ca. Dataset accessed 9/26/2017.

# Michigan shapefile
mich_shp <- shapefile("Data/GIS/Michigan_NoIsleRoyale.shp")

# LAGOS NE lakes (to account for border lakes, contains all lakes in Michigan or within 10km of Michigan)
lakes_4ha <- shapefile("C:/Ian_GIS/LAGOS-NE-GISv1.0/LAGOS_NE_All_Lakes_4ha/LAGOS_NE_4ha_within10km_mich.shp")

##### Main program ####
# Calculated as if lakes are patches in the same way as terrestrial patches
# not taking into account dispersal abilities of species
# essentially a landscape structure analysis based on lake shape and other patch-type stats
# not all metrics eventually used in developing lake connectivity indices

## convert lake subset polygon to raster
# Set up a raster "template" for a 30 m grid
rlagos <- raster(extent(lakes_4ha), res=30)

## Rasterize the shapefile (assigning NHD COMID as patch value)
mich_lake_raster <- rasterize(lakes_4ha , rlagos, field=as.numeric(lakes_4ha@data$lagoslakei))

# a simpler, vector-based alternative: https://rpubs.com/dgolicher/9458
# can also check out landscapemetrics package, or spatialEco::land.metrics
# calculate various patch metrics, PatchStat help file defines column names
# can join to shapefile by patchID
# PatchStat function seems great, but not clear how "core" is defined. Nothing obvious in function code
# source code has cryptic (to me) C code
# some back-calculations seem to be about 30m inward, but not sure if this is always 30m or based on cell width

test <- PatchStat(mat=mich_lake_raster, cellsize=30, latlon=F) 
test$edge.area <- test$area - test$core.area
test$edge.area.index <- test$edge.area/test$area

# Basic plots of PatchStat output
par(mfrow=c(3,3))
hist(test$area, main='Area', xlab='sq m')
hist(test$perimeter, main='Perimeter', xlab='m')
hist(test$perim.area.ratio, main='Perim/Area', xlab='m/sq m') #but decreases with increasing lake size
hist(test$core.area, main='Core area', xlab='sq m')
hist(test$edge.area, main='Edge area', xlab='sq m')
hist(test$shape.index, main='Shape index', xlab='Lo-Hi') #1=square shape, higher number=more complex shape. Corrects for perim/area issue
hist(test$frac.dim.index, main='Fractal dim index', xlab='Lo-Hi')#1=square, 2=highly complex shape. Also corrects for perim/area issue
hist(test$core.area.index, main='Core area index', xlab='Lo-Hi')#pct of patch that is core area
hist(test$edge.area.index, main='Edge area index', xlab='Lo-Hi')#pct of patch that is edge area

# save output
#write.csv(test, file="Data/MichiganLakePatchStats_wBorderStates.csv")

